#!/usr/bin/env python3
"""Find A223 neutrophils and score transferred Tumor HIF θ (persistent vs reverted)."""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
E27 = Path("/ix1/ylee/kor11/A223/E27/ocm/E27_gex_adt_ocm.h5ad")
E29 = Path("/ix1/ylee/kor11/A223/E29/ocm/E29_gex_adt_ocm.h5ad")
PALAK = Path("/ix1/ylee/Palak/A223/E27_E29_cell_type_metadata.csv")
GENE_MAP = Path("/ix1/ylee/Palak/MC38/index/grcm39_splici/index/gene_id_to_name.tsv")
OUT = Path("/ix1/ylee/kor11/A223/neutrophil_hypoxia")
THETA_QC = HERE / "tumor_theta_gene_qc.csv"

TAM_GENES = ("C1qa", "C1qb", "C1qc", "Adgre1", "Csf1r", "Cd68", "Fcgr1", "Trem2", "Mrc1", "Mertk")
NEU_GENES = ("S100a8", "S100a9", "Cxcr2", "Csf3r", "Retnlg", "Lcn2")
PALAK_SAMPLE = {
    "hypoxia_pos": "hypoxia_plus",
    "lactate_pos": "lactate_plus",
    "DN": "DN",
    "DP": "DP",
}
T_LOW, T_HIGH = 0.30, 0.70
NEU_MIN_UMI = 1500.0
TUMOR_MIN_UMI = 5000.0


def load_symbols() -> dict[str, str]:
    m = {}
    for line in GENE_MAP.read_text().splitlines():
        gid, name = line.split("\t", 1)
        m[gid] = name
    return m


def set_symbols(adata: ad.AnnData, mapping: dict[str, str]) -> None:
    names = [mapping.get(str(g), str(g)) for g in adata.var_names]
    adata.var["ensembl"] = np.asarray(adata.var_names.astype(str))
    adata.var_names = pd.Index(names)
    adata.var_names_make_unique()


def log_cp10k(X):
    X = X.astype(np.float64)
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    lib = np.asarray(X.sum(axis=1)).ravel()
    lib[lib <= 0] = 1.0
    X = X.multiply(1e4 / lib[:, None])
    X.data = np.log1p(X.data)
    return X.tocsr()


def gene_mean(X, var_names: np.ndarray, genes: tuple[str, ...]) -> np.ndarray:
    cols = [i for i, g in enumerate(var_names) if g in set(genes)]
    if not cols:
        return np.zeros(X.shape[0], dtype=np.float64)
    return np.asarray(X[:, cols].mean(axis=1)).ravel()


def csr_cols(mat, cols: list[int]) -> np.ndarray:
    if not cols:
        return np.zeros((mat.shape[0], 0), dtype=np.float64)
    sub = mat[:, cols]
    if sparse.issparse(sub):
        return np.asarray(sub.todense(), dtype=np.float64)
    return np.asarray(sub, dtype=np.float64)


def size_normalize(counts: np.ndarray, lib: np.ndarray) -> np.ndarray:
    lib = np.where(np.asarray(lib, dtype=np.float64).ravel() <= 0, 1.0, np.asarray(lib, dtype=np.float64).ravel())
    target = float(np.median(lib))
    return counts * (target / lib)[:, None]


def logistic_theta(h: np.ndarray, h_lo: float, h_hi: float) -> np.ndarray:
    span = max(h_hi - h_lo, 1e-6)
    h0 = 0.5 * (h_lo + h_hi)
    tau = span / 6.0
    return 1.0 / (1.0 + np.exp((h - h0) / tau))


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    if a.size < 5 or b.size < 5:
        return float("nan")
    sp = np.sqrt(((a.size - 1) * a.var(ddof=1) + (b.size - 1) * b.var(ddof=1)) / max(a.size + b.size - 2, 1))
    if sp < 1e-12:
        return 0.0
    return float((a.mean() - b.mean()) / sp)


def theta_state(theta: np.ndarray) -> np.ndarray:
    out = np.array(["partial"] * theta.size, dtype=object)
    out[theta <= T_LOW] = "persistent_like"
    out[theta >= T_HIGH] = "reverted_like"
    return out


def load_merged(mapping: dict[str, str]) -> ad.AnnData:
    parts = []
    for lane, path in (("E27", E27), ("E29", E29)):
        a = ad.read_h5ad(path)
        set_symbols(a, mapping)
        a.obs["lane"] = lane
        a.obs_names = pd.Index([f"{lane}_{b}" for b in a.obs_names.astype(str)])
        parts.append(a)
    adata = ad.concat(parts, join="inner", merge="same")
    adata.obs["barcode16"] = adata.obs["barcodes"].astype(str)
    return adata


def palak_table() -> pd.DataFrame:
    df = pd.read_csv(PALAK)
    parts = df["cell_barcode"].astype(str).str.extract(r"^(?P<bc>[ACGTN]+)-1-(?P<lane>e27|e29)_(?P<psample>.+)$")
    df = df.join(parts)
    df["lane"] = df["lane"].str.upper()
    df["sample_id"] = df["psample"].map(PALAK_SAMPLE)
    df["palak_neutrophil"] = df["cell_type_broad"].astype(str).str.strip().eq("Neutrophil")
    df["palak_tumor"] = df["cell_type_broad"].astype(str).str.strip().eq("Tumor")
    df["join"] = df["lane"] + "_" + df["bc"]
    return df


def annotate(adata: ad.AnnData) -> pd.DataFrame:
    var_names = np.asarray(adata.var_names.astype(str))
    X = log_cp10k(adata.layers["spliced"])
    lib_s = np.asarray(adata.layers["spliced"].sum(axis=1), dtype=np.float64).ravel()
    scores = pd.DataFrame(
        {
            "lane": adata.obs["lane"].astype(str).values,
            "sample_id": adata.obs["sample_id"].astype(str).values,
            "barcode16": adata.obs["barcode16"].astype(str).values,
            "gex_counts": np.asarray(adata.obs["gex_counts"], dtype=np.float64),
            "spliced_umi": lib_s,
            "ptprc": gene_mean(X, var_names, ("Ptprc",)),
            "tam": gene_mean(X, var_names, TAM_GENES),
            "neu": gene_mean(X, var_names, NEU_GENES),
            "s100a8": gene_mean(X, var_names, ("S100a8",)),
            "s100a9": gene_mean(X, var_names, ("S100a9",)),
            "cxcr2": gene_mean(X, var_names, ("Cxcr2",)),
        },
        index=adata.obs_names,
    )
    neu_cell = (scores["neu"] > scores["tam"] + 0.4) & (scores["neu"] > 1.5) & (scores["ptprc"] > 0.4)
    tumor_cell = (scores["ptprc"] < 0.25) & (scores["gex_counts"] >= 2000) & (scores["tam"] < 0.8) & (scores["neu"] < 1.0)
    scores["marker_group"] = "other"
    scores.loc[tumor_cell, "marker_group"] = "Tumor"
    scores.loc[neu_cell, "marker_group"] = "neutrophil"
    palak = palak_table()
    palak_idx = palak.drop_duplicates("join").set_index("join")
    key = scores["lane"] + "_" + scores["barcode16"]
    scores["palak_cell_type"] = key.map(palak_idx["cell_type"])
    scores["palak_broad"] = key.map(palak_idx["cell_type_broad"].astype(str).str.strip())
    scores["palak_neutrophil"] = key.map(palak_idx["palak_neutrophil"]).fillna(False).astype(bool)
    scores["palak_tumor"] = key.map(palak_idx["palak_tumor"]).fillna(False).astype(bool)
    scores["palak_matched"] = scores["palak_broad"].notna()
    scores["cell_group"] = scores["marker_group"]
    scores.loc[scores["palak_neutrophil"], "cell_group"] = "neutrophil"
    scores.loc[scores["palak_tumor"] & ~scores["palak_neutrophil"], "cell_group"] = "Tumor"
    return scores


def score_hif(adata: ad.AnnData, scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    qc = pd.read_csv(THETA_QC)
    use = qc["use_theta"].astype(str).str.lower().isin(["true", "1"])
    hif_genes = [g for g in qc.loc[use, "gene"].astype(str) if g in set(adata.var_names.astype(str))]
    w = qc.set_index("gene").loc[hif_genes, "cohens_d_e15_vs_e14"].clip(lower=0).to_numpy(dtype=np.float64)
    w = w / w.sum()
    var_pos = {str(g): i for i, g in enumerate(adata.var_names.astype(str))}
    cols = [var_pos[g] for g in hif_genes]
    lib_s = scores["spliced_umi"].to_numpy()
    spl = size_normalize(csr_cols(adata.layers["spliced"], cols), lib_s)
    log_s = np.log1p(spl)

    tumor_ref = (scores["cell_group"].eq("Tumor") & (lib_s >= TUMOR_MIN_UMI) & scores["sample_id"].eq("DN")).to_numpy()
    if int(tumor_ref.sum()) < 30:
        tumor_ref = (
            scores["cell_group"].eq("Tumor")
            & (lib_s >= TUMOR_MIN_UMI)
            & scores["sample_id"].isin(["DN", "lactate_plus"])
        ).to_numpy()
    tumor_hyp = (scores["cell_group"].eq("Tumor") & (lib_s >= TUMOR_MIN_UMI) & scores["sample_id"].eq("hypoxia_plus")).to_numpy()
    if tumor_ref.sum() < 15 or tumor_hyp.sum() < 15:
        raise SystemExit(f"tumor anchors too small: DN-like={tumor_ref.sum()} hypoxia+={tumor_hyp.sum()}")

    mu = log_s[tumor_ref].mean(axis=0)
    sd = log_s[tumor_ref].std(axis=0)
    sd = np.where(sd < 1e-12, 1.0, sd)
    z = (log_s - mu) / sd
    h = z @ w
    h_lo = float(np.median(h[tumor_ref]))
    h15q = float(np.quantile(h[tumor_hyp], 0.75))
    hi_mask = tumor_hyp & (h >= h15q)
    h_hi = float(np.median(h[hi_mask])) if hi_mask.any() else float(np.median(h[tumor_hyp]))
    theta = logistic_theta(h, h_lo, h_hi)
    state = theta_state(theta)

    out = scores.copy()
    out["hif_score"] = h
    out["theta_normoxic"] = theta
    out["hif_state"] = state
    out.loc[lib_s < NEU_MIN_UMI, "hif_state"] = "low_umi"

    gene_rows = []
    neu_ok = out["cell_group"].eq("neutrophil") & (lib_s >= NEU_MIN_UMI)
    for j, g in enumerate(hif_genes):
        hyp = log_s[neu_ok & out["sample_id"].eq("hypoxia_plus").to_numpy(), j]
        dn = log_s[neu_ok & out["sample_id"].eq("DN").to_numpy(), j]
        tum_h = log_s[tumor_hyp, j]
        tum_d = log_s[tumor_ref, j]
        gene_rows.append(
            {
                "gene": g,
                "weight": float(w[j]),
                "d_neu_hypoxia_vs_DN": cohens_d(hyp, dn),
                "d_tumor_hypoxia_vs_DN": cohens_d(tum_h, tum_d),
                "frac_nz": float((spl[:, j] > 0).mean()),
            }
        )
    gene_df = pd.DataFrame(gene_rows)
    meta = {
        "hif_genes": hif_genes,
        "h_lo_tumor_DN": h_lo,
        "h_hi_tumor_hypoxia": h_hi,
        "n_tumor_anchor_DN": int(tumor_ref.sum()),
        "n_tumor_anchor_hypoxia": int(tumor_hyp.sum()),
        "theta_low": T_LOW,
        "theta_high": T_HIGH,
        "neu_min_umi": NEU_MIN_UMI,
        "tumor_min_umi": TUMOR_MIN_UMI,
    }
    return out, gene_df, meta


def auroc_hif(cells: pd.DataFrame, group: str) -> dict:
    sub = cells[cells["cell_group"].eq(group) & cells["spliced_umi"].ge(NEU_MIN_UMI if group == "neutrophil" else TUMOR_MIN_UMI)]
    a = sub[sub["sample_id"].eq("hypoxia_plus")]
    b = sub[sub["sample_id"].eq("DN")]
    if len(a) < 8 or len(b) < 8:
        return {"n_hypoxia_plus": int(len(a)), "n_DN": int(len(b)), "auroc_hif": None}
    y = np.r_[np.ones(len(a)), np.zeros(len(b))]
    s = np.r_[a["hif_score"].to_numpy(), b["hif_score"].to_numpy()]
    return {
        "n_hypoxia_plus": int(len(a)),
        "n_DN": int(len(b)),
        "auroc_hif": float(roc_auc_score(y, s)),
        "median_theta_hypoxia_plus": float(a["theta_normoxic"].median()),
        "median_theta_DN": float(b["theta_normoxic"].median()),
        "median_theta_DP": float(sub.loc[sub["sample_id"].eq("DP"), "theta_normoxic"].median())
        if sub["sample_id"].eq("DP").any()
        else None,
    }


def plot_figures(cells: pd.DataFrame, outdir: Path) -> None:
    neu = cells[cells["cell_group"].eq("neutrophil") & cells["spliced_umi"].ge(NEU_MIN_UMI)]
    tum = cells[cells["cell_group"].eq("Tumor") & cells["spliced_umi"].ge(TUMOR_MIN_UMI)]
    order = ["hypoxia_plus", "DP", "DN", "lactate_plus"]
    states = ["persistent_like", "partial", "reverted_like"]
    colors = {"persistent_like": "#b2182b", "partial": "#f4a582", "reverted_like": "#2166ac"}

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    x = np.arange(len(order))
    bottom = np.zeros(len(order))
    for st in states:
        h = []
        for s in order:
            n = int((neu["sample_id"].eq(s)).sum())
            k = int(((neu["sample_id"].eq(s)) & (neu["hif_state"].eq(st))).sum())
            h.append(k / n if n else 0.0)
        ax.bar(x, h, bottom=bottom, color=colors[st], label=st.replace("_", " "))
        bottom += np.array(h)
    ax.set_xticks(x)
    ax.set_xticklabels(["hypoxia+", "DP", "DN", "lactate+"])
    ax.set_ylabel("fraction of neutrophils")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False)
    ax.set_title("Neutrophil HIF θ vs OCM stain")
    fig.tight_layout()
    fig.savefig(outdir / "neu_state_by_ocm.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    data, labs, cols = [], [], []
    for grp, df, c in (("neu", neu, "#1b9e77"), ("tumor", tum, "#d95f02")):
        for s in order:
            v = df.loc[df["sample_id"].eq(s), "theta_normoxic"].dropna().to_numpy()
            data.append(v)
            labs.append(f"{grp}\n{s.replace('_plus', '+')}")
            cols.append(c)
    bp = ax.boxplot(data, labels=labs, showfliers=False, patch_artist=True)
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c)
        patch.set_alpha(0.55)
    ax.axhline(T_LOW, color="#b2182b", ls="--", lw=0.8)
    ax.axhline(T_HIGH, color="#2166ac", ls="--", lw=0.8)
    ax.set_ylabel("θ (high = HIF-low / reverted-like)")
    ax.set_title("Transferred Tumor HIF θ")
    fig.tight_layout()
    fig.savefig(outdir / "theta_neu_vs_tumor.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    ax.scatter(neu["neu"], neu["tam"], c="#1b9e77", s=8, alpha=0.45, label="neutrophil")
    oth = cells[cells["cell_group"].eq("other")]
    ax.scatter(oth["neu"], oth["tam"], c="#bbbbbb", s=4, alpha=0.15, label="other")
    ax.set_xlabel("neutrophil marker (log CP10k)")
    ax.set_ylabel("TAM marker (log CP10k)")
    ax.legend(frameon=False, markerscale=2)
    fig.tight_layout()
    fig.savefig(outdir / "neu_vs_tam_markers.png", dpi=160)
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    mapping = load_symbols()
    adata = load_merged(mapping)
    scores = annotate(adata)
    cells, gene_df, meta = score_hif(adata, scores)

    neu = cells[cells["cell_group"].eq("neutrophil")]
    neu_qc = neu[neu["spliced_umi"] >= NEU_MIN_UMI]
    summary = {
        **meta,
        "n_cells": int(len(cells)),
        "n_palak_matched": int(cells["palak_matched"].sum()),
        "n_neutrophil_any": int(len(neu)),
        "n_neutrophil_qc": int(len(neu_qc)),
        "n_neutrophil_palak": int(cells["palak_neutrophil"].sum()),
        "n_neutrophil_marker": int(cells["marker_group"].eq("neutrophil").sum()),
        "neutrophil_by_ocm": neu_qc.groupby("sample_id").size().to_dict(),
        "neutrophil_state_by_ocm": neu_qc.groupby(["sample_id", "hif_state"]).size().unstack(fill_value=0).to_dict(),
        "tumor_state_by_ocm": cells[cells["cell_group"].eq("Tumor") & cells["spliced_umi"].ge(TUMOR_MIN_UMI)]
        .groupby(["sample_id", "hif_state"])
        .size()
        .unstack(fill_value=0)
        .to_dict(),
        "neutrophil_hif_auroc": auroc_hif(cells, "neutrophil"),
        "tumor_hif_auroc": auroc_hif(cells, "Tumor"),
    }
    cells.to_csv(OUT / "cells.csv")
    gene_df.to_csv(OUT / "hif_gene_qc.csv", index=False)
    neu_qc.to_csv(OUT / "neutrophils.csv")
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    plot_figures(cells, OUT)

    print("neutrophils (any / QC UMI>=1500):", len(neu), len(neu_qc))
    print(neu_qc.groupby(["sample_id", "hif_state"]).size().unstack(fill_value=0))
    print("neu HIF AUROC hypoxia+ vs DN:", summary["neutrophil_hif_auroc"])
    print("tumor HIF AUROC hypoxia+ vs DN:", summary["tumor_hif_auroc"])
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
