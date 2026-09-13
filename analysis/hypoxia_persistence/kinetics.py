#!/usr/bin/env python3
"""Tumor 1-D hypoxia → reoxygenation kinetics (spliced / unspliced).

Only `cell_group == Tumor`. E14S is the never-hypoxic anchor. E15S is a
mixture of persistent / reverting / reverted / memory. There is no GFP and
no time course; the derivative is identified from unspliced lag.

Pipeline
--------
1. Size-normalize spliced and unspliced by the **same** spliced UMI total.
   Separate u/s library scaling was destroying gene-level u/s (E14 tumor
   unspliced fraction is ~2× E15). Tumor cells below min_umi are dropped.
2. Tirosh S / G2M scores (mean log1p). Reported as covariates; they are
   **not** subtracted from HIF genes for θ (glycolytic HIF targets are
   collinear with growth, and residualizing them deletes the hypoxia axis).
3. Down-genes enter θ if Cohen's d of log1p size-normalized spliced
   (E15 vs E14) ≥ d_min. h = d-weighted mean of E14-z-scored log1p spliced.
4. θ = 1 / (1 + exp((h - h0) / τ)) with h0 midpoint of E14 median and
   E15-upper-quartile median, τ = |h_hyp - h_E14| / 6.
5. κ_g = E[u]/E[s] on **E15 persistent** cells only (lowest 20% θ). E14 is
   not mixed in: capture of unspliced differs by library. Ratio of means
   is used because unspliced is zero-inflated (median u/s = 0). Genes with
   κ_E14 / κ_persist above kappa_ratio_max stay in θ but drop from velocity.
6. Neighbor-smooth u and s on a HIF-spliced PCA graph, then
   v_g = u − κ s. Scale each gene by the persistent-cell MAD (same depth as
   E15), not E14 SD. v_reox = −weighted mean; high = transcription already
   down. Center so persist median is 0. Cycle residualization is E14-only
   and optional.
7. v_cut is the E15-persistent 95th percentile of |v_reox| (depth-matched
   SS null). Reverting requires hypoxic spliced (θ ≤ t_high) and inducing
   requires θ ≥ t_low, so high-unspliced outliers at the wrong end of θ
   are not called transients. E14 still gets the same |v| gate in the
   control table (unconstrained) as a capture/noise diagnostic.

Memory: Muc1 if detected, else Sod2, z vs E14. High θ + high memory
= memory, not full reversion.

HIF-1α protein is not observed; Hif1a mRNA is excluded from θ.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

HERE = Path(__file__).resolve().parent
GENESET_DIR = HERE / "genesets"
DEFAULT_H5AD = Path("/ix1/ylee/kor11/MC38/E14SE15S/E14SE15S_gex_adt_placed.h5ad")
DEFAULT_GENES = HERE / "reversion_key_genes.tsv"
NORMOXIA = "E14S"
HYPOXIA = "E15S"
TUMOR = "Tumor"
ALIASES = {"Ero1l": "Ero1a"}


def load_reversion_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["mouse"] = df["mouse"].map(lambda g: ALIASES.get(str(g), str(g)))
    df["sign"] = np.where(df["direction_on_reversion"].astype(str).str.startswith("down"), -1, 1)
    df.loc[df["mouse"].eq("Hif1a"), "sign"] = 0
    df.loc[df["mouse"].eq("Muc1"), "module"] = "memory"
    df.loc[df["mouse"].eq("Sod2"), "module"] = "memory"
    df.loc[df["module"].isna() & (df["sign"] < 0), "module"] = "hif_down"
    df.loc[df["module"].isna() & (df["sign"] > 0), "module"] = "reox_up"
    return df


def gene_index(var_names: np.ndarray, symbols: list[str]) -> dict[str, int]:
    pos = {str(g): i for i, g in enumerate(var_names)}
    return {s: pos[s] for s in symbols if s in pos}


def csr_cols(mat, cols: list[int]) -> np.ndarray:
    if not cols:
        return np.zeros((mat.shape[0], 0), dtype=np.float64)
    sub = mat[:, cols]
    if sparse.issparse(sub):
        return np.asarray(sub.todense(), dtype=np.float64)
    return np.asarray(sub, dtype=np.float64)


def size_normalize(counts: np.ndarray, lib: np.ndarray | None = None) -> np.ndarray:
    if lib is None:
        lib = counts.sum(axis=1)
    lib = np.where(np.asarray(lib, dtype=np.float64).ravel() <= 0, 1.0, np.asarray(lib, dtype=np.float64).ravel())
    target = float(np.median(lib))
    return counts * (target / lib)[:, None]


def knn_smooth(mat: np.ndarray, coords: np.ndarray, k: int) -> np.ndarray:
    n = mat.shape[0]
    k = int(max(2, min(k, n)))
    nn = NearestNeighbors(n_neighbors=k).fit(coords)
    ind = nn.kneighbors(return_distance=False)
    return mat[ind].mean(axis=1)


def mad(x: np.ndarray) -> float:
    med = float(np.median(x))
    return float(np.median(np.abs(x - med)) * 1.4826)


def mean_log1p_score(norm: np.ndarray, var_names: np.ndarray, symbols: list[str]) -> np.ndarray:
    idx = gene_index(var_names, symbols)
    if len(idx) < 5:
        raise SystemExit(f"cycle score: only {len(idx)} genes found")
    cols = list(idx.values())
    x = np.log1p(norm[:, cols])
    return x.mean(axis=1)


def ols_cycle_e14(y: np.ndarray, s: np.ndarray, g2m: np.ndarray, e14: np.ndarray, r2_min: float = 0.05):
    """Fit log1p ~ S + G2M on E14; return (b_s, b_g2m) or zeros if weak."""
    y14 = y[e14]
    X = np.column_stack([np.ones(e14.sum()), s[e14], g2m[e14]])
    coef, *_ = np.linalg.lstsq(X, y14, rcond=None)
    pred = X @ coef
    sst = np.sum((y14 - y14.mean()) ** 2)
    r2 = 1.0 - np.sum((y14 - pred) ** 2) / sst if sst > 0 else 0.0
    if r2 < r2_min:
        return 0.0, 0.0, float(r2)
    return float(coef[1]), float(coef[2]), float(r2)


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    n1, n2 = a.size, b.size
    if n1 < 5 or n2 < 5:
        return 0.0
    v1, v2 = a.var(ddof=1), b.var(ddof=1)
    sp = np.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / max(n1 + n2 - 2, 1))
    if sp < 1e-12:
        return 0.0
    return float((a.mean() - b.mean()) / sp)


def logistic_theta(h: np.ndarray, h_lo: float, h_hi: float) -> np.ndarray:
    """Map score h (high = hypoxic) to θ in (0,1); anchors → ~0.05 and ~0.95."""
    span = max(h_hi - h_lo, 1e-6)
    h0 = 0.5 * (h_lo + h_hi)
    tau = span / 6.0
    return 1.0 / (1.0 + np.exp((h - h0) / tau))


def z_vs_ref(x: np.ndarray, ref_mask: np.ndarray) -> np.ndarray:
    mu = x[ref_mask].mean()
    sd = x[ref_mask].std()
    sd = sd if sd > 1e-12 else 1.0
    return (x - mu) / sd


def classify(
    sample: np.ndarray,
    theta: np.ndarray,
    v_reox: np.ndarray,
    memory_z: np.ndarray,
    v_cut: float,
    t_low: float,
    t_high: float,
    mem_cut: float,
) -> np.ndarray:
    """Velocity only counts if θ still matches that direction of travel."""
    n = sample.size
    out = np.array(["unassigned"] * n, dtype=object)
    e14 = sample == NORMOXIA
    e15 = sample == HYPOXIA
    toward_normoxia = (v_reox >= v_cut) & (theta <= t_high)
    toward_hypoxia = (v_reox <= -v_cut) & (theta >= t_low)
    out[toward_normoxia] = "reverting"
    out[toward_hypoxia] = "inducing"
    moving = toward_normoxia | toward_hypoxia
    out[e14 & ~moving] = "never_hypoxic"
    rest = e15 & ~moving
    mem = rest & (theta >= t_high) & (memory_z >= mem_cut)
    out[mem] = "memory"
    rest = rest & ~mem
    out[rest & (theta >= t_high)] = "reverted"
    out[rest & (theta <= t_low)] = "persistent"
    out[rest & (theta > t_low) & (theta < t_high)] = "partial"
    return out


def direction_control(sample: np.ndarray, v_reox: np.ndarray, v_cut: float, theta: np.ndarray, t_low: float, t_high: float) -> pd.DataFrame:
    rows = []
    for name, mask in (("E14S", sample == NORMOXIA), ("E15S", sample == HYPOXIA)):
        n = int(mask.sum())
        n_rev = int((mask & (v_reox >= v_cut)).sum())
        n_ind = int((mask & (v_reox <= -v_cut)).sum())
        n_rev_ph = int((mask & (v_reox >= v_cut) & (theta <= t_high)).sum())
        n_ind_ph = int((mask & (v_reox <= -v_cut) & (theta >= t_low)).sum())
        rows.append(
            {
                "sample": name,
                "n": n,
                "v_cut": v_cut,
                "n_reverting": n_rev_ph,
                "n_inducing": n_ind_ph,
                "n_reverting_unconstrained": n_rev,
                "n_inducing_unconstrained": n_ind,
                "frac_reverting": n_rev_ph / n if n else np.nan,
                "frac_inducing": n_ind_ph / n if n else np.nan,
                "median_v_reox": float(np.median(v_reox[mask])),
                "q05_v_reox": float(np.quantile(v_reox[mask], 0.05)),
                "q95_v_reox": float(np.quantile(v_reox[mask], 0.95)),
            }
        )
    return pd.DataFrame(rows)


def run_kinetics(
    adata,
    genes_path: Path,
    d_min: float = 0.20,
    min_spl_nz: float = 0.15,
    min_uns_nz: float = 0.10,
    min_umi: float = 5000,
    kappa_ratio_max: float = 4.0,
    t_low: float = 0.30,
    t_high: float = 0.70,
    mem_cut: float = 1.0,
    residualize_v_cycle: bool = True,
    smooth_k: int = 30,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if "cell_group" not in adata.obs:
        raise SystemExit("obs['cell_group'] missing; run annotate_cell_groups.py first")
    table = load_reversion_table(genes_path)
    var_names = np.asarray(adata.var_names.astype(str))
    tumor = adata.obs["cell_group"].astype(str).to_numpy() == TUMOR
    if tumor.sum() < 50:
        raise SystemExit(f"too few Tumor cells: {tumor.sum()}")

    ad_t = adata[tumor].copy()
    sample = ad_t.obs["sample"].astype(str).to_numpy()
    e14 = sample == NORMOXIA
    e15 = sample == HYPOXIA
    if e14.sum() < 20 or e15.sum() < 20:
        raise SystemExit(f"Tumor needs both samples; E14={e14.sum()} E15={e15.sum()}")

    lib_s = np.asarray(ad_t.layers["spliced"].sum(axis=1), dtype=np.float64).ravel()
    lib_u = np.asarray(ad_t.layers["unspliced"].sum(axis=1), dtype=np.float64).ravel()
    keep = lib_s >= min_umi
    if keep.sum() < 50 or (keep & e14).sum() < 15:
        raise SystemExit(
            f"min_umi={min_umi} leaves Tumor n={keep.sum()} E14={(keep & e14).sum()}; "
            "E14 tumor libraries are much smaller than E15"
        )
    dropped = int((~keep).sum())
    ad_t = ad_t[keep].copy()
    sample = sample[keep]
    e14 = e14[keep]
    e15 = e15[keep]
    lib_s = lib_s[keep]
    lib_u = lib_u[keep]
    print(f"Tumor QC min spliced UMI {min_umi}: kept {int(keep.sum())} dropped {dropped} (E14={e14.sum()} E15={e15.sum()})")
    print(
        f"median spliced UMI E14={np.median(lib_s[e14]):.0f} E15={np.median(lib_s[e15]):.0f}; "
        f"median U/S E14={np.median(lib_u[e14] / np.maximum(lib_s[e14], 1)):.3f} "
        f"E15={np.median(lib_u[e15] / np.maximum(lib_s[e15], 1)):.3f}"
    )
    cycle_spec = json.loads((GENESET_DIR / "tirosh_cell_cycle_mouse.json").read_text())
    cycle_syms = cycle_spec["S"] + cycle_spec["G2M"]
    cix = gene_index(var_names, cycle_syms)
    cycle_norm = size_normalize(csr_cols(ad_t.layers["spliced"], [cix[g] for g in cix]), lib_s)
    cycle_names = np.array(list(cix))
    s_score = mean_log1p_score(cycle_norm, cycle_names, cycle_spec["S"])
    g2m_score = mean_log1p_score(cycle_norm, cycle_names, cycle_spec["G2M"])
    s0, g0 = float(s_score[e14].mean()), float(g2m_score[e14].mean())

    symbols = table["mouse"].tolist()
    gix = gene_index(var_names, symbols)
    spl = size_normalize(csr_cols(ad_t.layers["spliced"], [gix[g] for g in gix]), lib_s)
    uns = size_normalize(csr_cols(ad_t.layers["unspliced"], [gix[g] for g in gix]), lib_s)
    genes = list(gix)
    col = {g: i for i, g in enumerate(genes)}

    rows = []
    for g in genes:
        j = col[g]
        meta = table.loc[table["mouse"] == g].iloc[0]
        log_s = np.log1p(spl[:, j])
        _bs, _bg, r2s = ols_cycle_e14(log_s, s_score, g2m_score, e14)
        _bu, _bgu, r2u = ols_cycle_e14(np.log1p(uns[:, j]), s_score, g2m_score, e14)
        d = cohens_d(log_s[e15], log_s[e14])
        rows.append(
            {
                "gene": g,
                "module": meta["module"],
                "rank": int(meta["rank"]),
                "spliced_frac": float((spl[:, j] > 0).mean()),
                "unspliced_frac": float((uns[:, j] > 0).mean()),
                "cohens_d_e15_vs_e14": d,
                "cycle_r2_spliced_e14": r2s,
                "cycle_r2_unspliced_e14": r2u,
            }
        )
    qc = pd.DataFrame(rows)

    hif = qc["module"].eq("hif_down") & (qc["spliced_frac"] >= min_spl_nz) & (qc["cohens_d_e15_vs_e14"] >= d_min)
    if int(hif.sum()) < 4:
        raise SystemExit(f"only {int(hif.sum())} HIF-down genes passed QC")
    hif_genes = qc.loc[hif, "gene"].tolist()
    weights = qc.set_index("gene").loc[hif_genes, "cohens_d_e15_vs_e14"].to_numpy(dtype=np.float64)
    weights = weights / weights.sum()

    z_s = np.column_stack([z_vs_ref(np.log1p(spl[:, col[g]]), e14) for g in hif_genes])
    h = z_s @ weights
    h14 = float(np.median(h[e14]))
    h15q = float(np.quantile(h[e15], 0.75))
    h_hyp = float(np.median(h[e15][h[e15] >= h15q])) if np.any(h[e15] >= h15q) else float(np.median(h[e15]))
    theta = logistic_theta(h, h14, h_hyp)

    persist_ss = e15 & (theta <= np.quantile(theta[e15], 0.20))
    if persist_ss.sum() < 20:
        raise SystemExit("too few E15 persistent cells to estimate κ")

    def kappa_means(j: int, mask: np.ndarray) -> float:
        if mask.sum() < 10:
            return np.nan
        sm = float(spl[mask, j].mean())
        if sm < 1e-6:
            return np.nan
        return float(uns[mask, j].mean() / (sm + 1e-6))

    kappa = {}
    for g in hif_genes:
        j = col[g]
        k14, k15 = kappa_means(j, e14), kappa_means(j, persist_ss)
        ratio = np.nan
        if np.isfinite(k14) and np.isfinite(k15) and min(k14, k15) > 1e-8:
            ratio = max(k14, k15) / min(k14, k15)
        qc.loc[qc["gene"] == g, "kappa_e14"] = k14
        qc.loc[qc["gene"] == g, "kappa_persist"] = k15
        qc.loc[qc["gene"] == g, "kappa_ss"] = k15
        qc.loc[qc["gene"] == g, "kappa_ratio"] = ratio
        uns_ok = float(qc.loc[qc["gene"] == g, "unspliced_frac"].iloc[0]) >= min_uns_nz
        ratio_ok = (not np.isfinite(ratio)) or (ratio <= kappa_ratio_max)
        use_v = uns_ok and np.isfinite(k15) and k15 > 1e-8 and ratio_ok
        qc.loc[qc["gene"] == g, "use_velocity"] = use_v
        if use_v:
            kappa[g] = k15

    qc["use_theta"] = qc["gene"].isin(hif_genes)
    qc["use_velocity"] = qc["use_velocity"].fillna(False)

    v_genes = list(kappa)
    if not v_genes:
        raise SystemExit("no genes passed unspliced/κ QC for velocity")
    spl_v = np.column_stack([spl[:, col[g]] for g in v_genes])
    uns_v = np.column_stack([uns[:, col[g]] for g in v_genes])
    k_vec = np.asarray([kappa[g] for g in v_genes], dtype=np.float64)
    if smooth_k and smooth_k > 1:
        pcs = PCA(n_components=min(8, max(2, len(hif_genes))), random_state=0).fit_transform(z_s)
        spl_v = knn_smooth(spl_v, pcs, smooth_k)
        uns_v = knn_smooth(uns_v, pcs, smooth_k)
        print(f"velocity moments: kNN k={smooth_k} on HIF spliced PCA")
    v_mat = uns_v - spl_v * k_vec
    gene_scale = np.array([mad(v_mat[persist_ss, j]) for j in range(len(v_genes))], dtype=np.float64)
    gene_scale = np.where(gene_scale < 1e-8, 1.0, gene_scale)
    v_w = qc.set_index("gene").loc[v_genes, "cohens_d_e15_vs_e14"].to_numpy(dtype=np.float64)
    v_w = v_w / v_w.sum()
    v_reox = -((v_mat / gene_scale) @ v_w)
    v_reox = v_reox - float(np.median(v_reox[persist_ss]))
    n_genes_down = np.sum(v_mat < 0, axis=1).astype(np.int32)
    cycle_coef = (0.0, 0.0)
    if residualize_v_cycle:
        X14 = np.column_stack([s_score[e14] - s0, g2m_score[e14] - g0])
        coef, *_ = np.linalg.lstsq(X14, v_reox[e14], rcond=None)
        cycle_coef = (float(coef[0]), float(coef[1]))
        v_reox = v_reox - coef[0] * (s_score - s0) - coef[1] * (g2m_score - g0)
        v_reox = v_reox - float(np.median(v_reox[persist_ss]))
    print(f"v_reox cycle residualization: {residualize_v_cycle} b_S={cycle_coef[0]:.3f} b_G2M={cycle_coef[1]:.3f}")

    mem_genes = []
    for g in ["Muc1", "Sod2"]:
        if g not in col:
            continue
        nz = float(qc.loc[qc["gene"] == g, "spliced_frac"].iloc[0]) if g in set(qc["gene"]) else 0.0
        if nz >= 0.05:
            mem_genes.append(g)
        qc.loc[qc["gene"] == g, "use_memory"] = nz >= 0.05
    if not mem_genes:
        mem_genes = [g for g in ["Sod2"] if g in col]
    mem_z = (
        np.mean(np.column_stack([z_vs_ref(np.log1p(spl[:, col[g]]), e14) for g in mem_genes]), axis=1)
        if mem_genes
        else np.zeros(e14.size)
    )

    v_cut = float(np.quantile(np.abs(v_reox[persist_ss]), 0.95))
    state = classify(sample, theta, v_reox, mem_z, v_cut, t_low, t_high, mem_cut)
    ctrl = direction_control(sample, v_reox, v_cut, theta, t_low, t_high)
    e14_ind = float(((sample == NORMOXIA) & (state == "inducing")).mean()) if e14.any() else 1.0
    e15_ind = float(((sample == HYPOXIA) & (state == "inducing")).mean()) if e15.any() else 0.0
    e14_rev = float(((sample == NORMOXIA) & (state == "reverting")).mean()) if e14.any() else 1.0
    e15_rev = float(((sample == HYPOXIA) & (state == "reverting")).mean()) if e15.any() else 0.0
    print(
        f"E15 vs E14 floors: reverting {e15_rev:.4f} vs {e14_rev:.4f}; "
        f"inducing {e15_ind:.4f} vs {e14_ind:.4f} (E15 labels kept even if not enriched; "
        "compare control table)"
    )

    cells = pd.DataFrame(
        {
            "sample": sample,
            "cell_group": TUMOR,
            "cycle_s": s_score,
            "cycle_g2m": g2m_score,
            "hif_score": h,
            "theta_normoxic": theta,
            "v_reox": v_reox,
            "v_hypoxia": -v_reox,
            "n_genes_down": n_genes_down,
            "memory_z": mem_z,
            "hypoxia_kinetics_state": state,
        },
        index=ad_t.obs_names,
    )
    if "barcodes" in ad_t.obs:
        cells["barcodes"] = ad_t.obs["barcodes"].astype(str).values
    if "X_umap" in ad_t.obsm:
        cells["umap_x"] = np.asarray(ad_t.obsm["X_umap"])[:, 0]
        cells["umap_y"] = np.asarray(ad_t.obsm["X_umap"])[:, 1]

    print(f"Tumor n={len(cells)} E14={e14.sum()} E15={e15.sum()}")
    print("HIF-down used for θ:", ", ".join(hif_genes))
    print("velocity genes:", ", ".join(kappa.keys()))
    print("memory genes:", ", ".join(mem_genes))
    print(
        f"θ anchors E14 median={h14:.3f} E15 persistent={h_hyp:.3f}  "
        f"|v| cut (E15 persist 95%)={v_cut:.3f}  persist n={int(persist_ss.sum())}"
    )
    print(cells.groupby(["sample", "hypoxia_kinetics_state"]).size().unstack(fill_value=0))
    print(cells.groupby("hypoxia_kinetics_state")[["theta_normoxic", "v_reox", "v_hypoxia", "memory_z", "cycle_s"]].median())
    print("direction control (v_cut from E15 persist; phenotype-constrained counts):")
    print(ctrl.to_string(index=False))
    return cells, qc, ctrl


def attach_to_full(adata, cells: pd.DataFrame) -> None:
    n = adata.n_obs
    idx = pd.Index(adata.obs_names)
    for col in ["cycle_s", "cycle_g2m", "hif_score", "theta_normoxic", "v_reox", "v_hypoxia", "memory_z"]:
        v = np.full(n, np.nan, dtype=np.float64)
        loc = idx.get_indexer(cells.index)
        v[loc] = cells[col].to_numpy(dtype=np.float64)
        adata.obs[col] = v
    if "n_genes_down" in cells.columns:
        nd = np.full(n, np.nan, dtype=np.float64)
        loc = idx.get_indexer(cells.index)
        nd[loc] = cells["n_genes_down"].to_numpy(dtype=np.float64)
        adata.obs["n_genes_down"] = nd
    st = np.array(["not_tumor"] * n, dtype=object)
    tumor = adata.obs["cell_group"].astype(str).to_numpy() == TUMOR
    st[tumor] = "tumor_low_umi"
    st[idx.get_indexer(cells.index)] = cells["hypoxia_kinetics_state"].to_numpy()
    adata.obs["hypoxia_kinetics_state"] = st


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", default=str(DEFAULT_H5AD))
    ap.add_argument("--genes", default=str(DEFAULT_GENES))
    ap.add_argument("--out-cells", default=str(HERE / "hypoxia_kinetics.csv"))
    ap.add_argument("--out-qc", default=str(HERE / "hypoxia_kinetics_gene_qc.csv"))
    ap.add_argument("--out-control", default=str(HERE / "hypoxia_kinetics_control.csv"))
    ap.add_argument("--write-h5ad", action="store_true")
    ap.add_argument("--smooth-k", type=int, default=30)
    args = ap.parse_args()

    import anndata as ad

    adata = ad.read_h5ad(args.h5ad)
    cells, qc, ctrl = run_kinetics(adata, Path(args.genes), smooth_k=args.smooth_k)
    Path(args.out_cells).parent.mkdir(parents=True, exist_ok=True)
    cells.to_csv(args.out_cells)
    qc.to_csv(args.out_qc, index=False)
    ctrl.to_csv(args.out_control, index=False)
    print(f"wrote {args.out_cells}, {args.out_qc}, {args.out_control}")
    if args.write_h5ad:
        attach_to_full(adata, cells)
        adata.write_h5ad(args.h5ad)
        print(f"updated {args.h5ad}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
