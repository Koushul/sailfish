#!/usr/bin/env python3
"""NADPH oxidase / priming / NRF2 RNA in DCF+ neutrophils (burst is post-transcriptional)."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
NEU_PY = HERE / "run.py"
E14E15 = Path("/ix1/ylee/kor11/MC38/E14SE15S/E14SE15S_gex_adt_placed.h5ad")
OUT = Path("/ix1/ylee/kor11/A223/neutrophil_hypoxia/ros_rna")
REPO = HERE
NEU_UMI = 1500.0

NOX = ("Cybb", "Cyba", "Ncf1", "Ncf2", "Ncf4", "Rac2")
PRIMING = ("Cxcl1", "Cxcl2", "Il1a", "Il1b", "Il1rn", "Icam1")
PRIMING_TNF = ("Ccl3", "Ccl4", "Tnf", "Nfkbia", "Tnfaip3")
NRF2 = ("Hmox1", "Nqo1", "Gclc", "Gclm", "Txnrd1", "Sod2", "Gpx1")
HIF = ("Ldha", "Vegfa", "Pgk1", "Eno1", "Slc2a1", "P4ha1")
MATURE = ("Fpr1", "Cxcr2", "C5ar1")
SETS = {
    "nox_capacity": NOX,
    "priming_wright2013": PRIMING,
    "priming_tnf": PRIMING_TNF,
    "nrf2_ros_response": NRF2,
    "hif_tumor_module": HIF,
    "mature_receptors": MATURE,
}


def load_neu():
    spec = importlib.util.spec_from_file_location("a223_neu_hif", NEU_PY)
    neu = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(neu)
    return neu


def module_score(S, lib, pos, genes):
    ix = [pos[g] for g in genes if g in pos]
    used = [g for g in genes if g in pos]
    if not ix:
        return np.zeros(S.shape[0]), used
    x = np.asarray(S[:, ix].sum(axis=1)).ravel()
    return np.log1p(x * (1e4 / np.maximum(lib, 1.0))), used


def gene_log1p(S, lib, j):
    raw = np.asarray(S[:, [j]].todense()).ravel()
    return np.log1p(raw * (1e4 / np.maximum(lib, 1.0)))


def contrast(score, pos_mask, neg_mask) -> dict:
    a, b = score[pos_mask], score[neg_mask]
    out = {
        "n_pos": int(a.size),
        "n_neg": int(b.size),
        "median_pos": float(np.median(a)) if a.size else np.nan,
        "median_neg": float(np.median(b)) if b.size else np.nan,
        "cohens_d": np.nan,
        "mwu_p": np.nan,
        "auroc": np.nan,
    }
    if min(a.size, b.size) < 8:
        return out
    sp = np.sqrt(0.5 * (a.var(ddof=1) + b.var(ddof=1)))
    out["cohens_d"] = float((a.mean() - b.mean()) / sp) if sp > 1e-12 else 0.0
    out["mwu_p"] = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    y = np.r_[np.ones(a.size), np.zeros(b.size)]
    out["auroc"] = float(roc_auc_score(y, np.r_[a, b]))
    return out


def plot_modules(df: pd.DataFrame, path: Path) -> None:
    sub = df[df["split"] == "A223 DCF+ vs DP neutrophils"].copy()
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    y = np.arange(len(sub))
    ax.barh(y, sub["auroc"].to_numpy() - 0.5, left=0.5, color="#2c3e50")
    ax.axvline(0.5, color="#7f8c8d", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["module"].tolist())
    ax.set_xlim(0.2, 0.8)
    ax.set_xlabel("AUROC (DCF+ > DP)")
    ax.set_title("A223 neutrophil RNA: DCF+ vs DP")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    neu = load_neu()
    mapping = neu.load_symbols()
    adata = neu.load_merged(mapping)
    scores = neu.annotate(adata)
    S = adata.layers["spliced"]
    if not hasattr(S, "tocsr"):
        from scipy import sparse

        S = sparse.csr_matrix(S)
    else:
        S = S.tocsr()
    lib = scores["spliced_umi"].to_numpy()
    pos = {str(g): i for i, g in enumerate(adata.var_names.astype(str))}
    neu_m = scores["palak_neutrophil"].to_numpy() & (lib >= NEU_UMI)
    gate = scores["sample_id"].astype(str).to_numpy()
    dcf = neu_m & (gate == "hypoxia_plus")
    dp = neu_m & (gate == "DP")
    dn = neu_m & (gate == "DN")

    rows, gene_rows = [], []
    scores_mod = {}
    for name, genes in SETS.items():
        sc, used = module_score(S, lib, pos, genes)
        scores_mod[name] = sc
        for split, pmask, nmask in (
            ("A223 DCF+ vs DN neutrophils", dcf, dn),
            ("A223 DCF+ vs DP neutrophils", dcf, dp),
            ("A223 DP vs DN neutrophils", dp, dn),
        ):
            met = contrast(sc, pmask, nmask)
            met.update({"split": split, "module": name, "genes": ",".join(used)})
            rows.append(met)
        for g in used:
            gl = gene_log1p(S, lib, pos[g])
            gene_rows.append(
                {
                    "module": name,
                    "gene": g,
                    "median_DCF": float(np.median(gl[dcf])),
                    "median_DP": float(np.median(gl[dp])),
                    "median_DN": float(np.median(gl[dn])) if dn.any() else np.nan,
                    "nz_DCF": float((np.asarray(S[dcf, pos[g]].todense()).ravel() > 0).mean()) if dcf.any() else np.nan,
                    "nz_DP": float((np.asarray(S[dp, pos[g]].todense()).ravel() > 0).mean()) if dp.any() else np.nan,
                    **{f"auroc_{k}": contrast(gl, a, b).get("auroc", np.nan) for k, a, b in (
                        ("DCF_vs_DP", dcf, dp),
                        ("DCF_vs_DN", dcf, dn),
                    )},
                }
            )

    ref = ad.read_h5ad(E14E15)
    neu.set_symbols(ref, mapping)
    Sr = ref.layers["spliced"]
    Sr = Sr.tocsr() if hasattr(Sr, "tocsr") else Sr
    libr = np.asarray(Sr.sum(axis=1), dtype=np.float64).ravel()
    posr = {str(g): i for i, g in enumerate(ref.var_names.astype(str))}
    sample = ref.obs["sample"].astype(str).to_numpy()
    cg = ref.obs["cell_group"].astype(str).to_numpy()
    e14n = (cg == "neutrophil") & (libr >= NEU_UMI) & (sample == "E14S")
    e15n = (cg == "neutrophil") & (libr >= NEU_UMI) & (sample == "E15S")
    for name, genes in SETS.items():
        sc, used = module_score(Sr, libr, posr, genes)
        met = contrast(sc, e15n, e14n)
        met.update({"split": "E15 vs E14 neutrophils", "module": name, "genes": ",".join(used)})
        rows.append(met)

    table = pd.DataFrame(rows)
    genes_tbl = pd.DataFrame(gene_rows)
    table.to_csv(OUT / "module_contrasts.csv", index=False)
    genes_tbl.to_csv(OUT / "gene_medians.csv", index=False)
    table.to_csv(REPO / "ros_rna_module_contrasts.csv", index=False)
    genes_tbl.to_csv(REPO / "ros_rna_gene_medians.csv", index=False)
    plot_modules(table, OUT / "dcf_vs_dp_module_auroc.png")
    show = table[["split", "module", "n_pos", "n_neg", "auroc", "cohens_d", "mwu_p", "median_pos", "median_neg"]]
    print(show.to_string(index=False))
    print("QC neutrophils", int(neu_m.sum()), "DCF", int(dcf.sum()), "DP", int(dp.sum()), "DN", int(dn.sum()))
    del adata, ref


if __name__ == "__main__":
    main()
