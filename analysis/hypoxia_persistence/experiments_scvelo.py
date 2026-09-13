#!/usr/bin/env python3
"""Compare scVelo (all cells) vs Tumor 1-D kinetics, with/without cell cycle."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from kinetics import (  # noqa: E402
    DEFAULT_GENES,
    DEFAULT_H5AD,
    GENESET_DIR,
    csr_cols,
    gene_index,
    mean_log1p_score,
    run_kinetics,
    size_normalize,
)

OUT = HERE / "experiments"
HIF_THETA = [
    "Car9", "Vegfa", "Ldha", "Pdk1", "Egln1", "Pgk1", "Eno1", "Aldoa",
    "Angptl4", "P4ha1", "Ero1a", "Ankrd37", "Serpine1", "Cited2",
]


def cycle_scores(adata) -> tuple[np.ndarray, np.ndarray]:
    spec = json.loads((GENESET_DIR / "tirosh_cell_cycle_mouse.json").read_text())
    vn = np.asarray(adata.var_names.astype(str))
    lib = np.asarray(adata.layers["spliced"].sum(axis=1), dtype=np.float64).ravel()
    cix = gene_index(vn, spec["S"] + spec["G2M"])
    mat = size_normalize(csr_cols(adata.layers["spliced"], [cix[g] for g in cix]), lib)
    names = np.array(list(cix))
    return mean_log1p_score(mat, names, spec["S"]), mean_log1p_score(mat, names, spec["G2M"])


def hif_expr_score(adata) -> np.ndarray:
    vn = np.asarray(adata.var_names.astype(str))
    lib = np.asarray(adata.layers["spliced"].sum(axis=1), dtype=np.float64).ravel()
    gix = gene_index(vn, HIF_THETA)
    mat = size_normalize(csr_cols(adata.layers["spliced"], [gix[g] for g in gix]), lib)
    return np.log1p(mat).mean(axis=1)


def residualize_counts(Y: np.ndarray, s: np.ndarray, g2m: np.ndarray, fit: np.ndarray) -> np.ndarray:
    X = np.column_stack([np.ones(s.size), s, g2m])
    Xf = X[fit]
    Ylog = np.log1p(np.asarray(Y, dtype=np.float64))
    out = np.empty_like(Ylog)
    for j in range(Ylog.shape[1]):
        coef, *_ = np.linalg.lstsq(Xf, Ylog[fit, j], rcond=None)
        cov = coef[1] * (s - s[fit].mean()) + coef[2] * (g2m - g2m[fit].mean())
        out[:, j] = np.expm1(np.clip(Ylog[:, j] - cov, 0, None))
    return out


def gene_set_velocity(adata, symbols: list[str]) -> np.ndarray:
    names = [str(g) for g in adata.var_names]
    pos = {g: i for i, g in enumerate(names)}
    V = adata.layers["velocity"]
    if sparse.issparse(V):
        V = V.toarray()
    cols = [pos[g] for g in symbols if g in pos]
    if not cols:
        return np.full(adata.n_obs, np.nan)
    block = np.asarray(V[:, cols], dtype=np.float64)
    return np.nanmean(block, axis=1)


def spearman(x, y) -> tuple[float, float]:
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 30:
        return np.nan, np.nan
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p)


def run_scvelo(adata, cycle: str, e14_mask: np.ndarray):
    import scvelo as scv

    scv.settings.verbosity = 1
    b = adata.copy()
    if "spliced" not in b.layers:
        raise SystemExit("need spliced/unspliced layers")
    spec = json.loads((GENESET_DIR / "tirosh_cell_cycle_mouse.json").read_text())
    scv.pp.filter_genes(b, min_shared_counts=20)
    scv.pp.normalize_per_cell(b)
    S = b.layers["spliced"]
    if sparse.issparse(S):
        S = S.toarray()
    logS = np.log1p(np.asarray(S, dtype=np.float64))
    logS[~np.isfinite(logS)] = 0
    var = logS.var(axis=0)
    top = set(np.argsort(var)[-2000:].tolist())
    names = [str(g) for g in b.var_names]
    pos = {g: i for i, g in enumerate(names)}
    for g in HIF_THETA + spec["S"] + spec["G2M"]:
        if g in pos:
            top.add(pos[g])
    b = b[:, sorted(top)].copy()
    if cycle != "none":
        fit = e14_mask if cycle == "e14" else np.ones(b.n_obs, dtype=bool)
        s = b.obs["cycle_s"].to_numpy(dtype=np.float64)
        g2m = b.obs["cycle_g2m"].to_numpy(dtype=np.float64)
        for layer in ("spliced", "unspliced"):
            Y = b.layers[layer]
            if sparse.issparse(Y):
                Y = Y.toarray()
            b.layers[layer] = residualize_counts(Y, s, g2m, fit)
            if layer == "spliced":
                b.X = b.layers[layer]
    scv.pp.moments(b, n_pcs=30, n_neighbors=30)
    scv.tl.velocity(b, mode="stochastic")
    return b


def summarize(run: str, v_hif, v_cyc, hif, s_score, v_reox=None) -> dict:
    r_hs, p_hs = spearman(v_hif, s_score)
    r_cs, p_cs = spearman(v_cyc, s_score)
    r_hh, p_hh = spearman(v_hif, hif)
    r_ours, p_ours = (np.nan, np.nan)
    if v_reox is not None:
        r_ours, p_ours = spearman(-v_hif, v_reox)
    mag_c = np.nanmean(np.abs(v_cyc))
    mag_h = np.nanmean(np.abs(v_hif))
    return {
        "run": run,
        "n": int(np.isfinite(v_hif).sum()),
        "spearman_vHIF_vs_S": r_hs,
        "p_vHIF_vs_S": p_hs,
        "spearman_vCycle_vs_S": r_cs,
        "p_vCycle_vs_S": p_cs,
        "spearman_vHIF_vs_hifExpr": r_hh,
        "p_vHIF_vs_hifExpr": p_hh,
        "spearman_minus_vHIF_vs_our_v_reox": r_ours,
        "p_minus_vHIF_vs_our_v_reox": p_ours,
        "mean_abs_vCycle": mag_c,
        "mean_abs_vHIF": mag_h,
        "cycle_over_hif_ratio": mag_c / mag_h if mag_h and np.isfinite(mag_h) and mag_h > 0 else np.nan,
        "frac_cycle_dominates": float(np.mean(np.abs(v_cyc) > np.abs(v_hif))),
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(DEFAULT_H5AD)
    s_score, g2m_score = cycle_scores(adata)
    adata.obs["cycle_s"] = s_score
    adata.obs["cycle_g2m"] = g2m_score
    adata.obs["hif_expr"] = hif_expr_score(adata)
    e14 = adata.obs["sample"].astype(str).to_numpy() == "E14S"
    spec = json.loads((GENESET_DIR / "tirosh_cell_cycle_mouse.json").read_text())
    cycle_genes = spec["S"] + spec["G2M"]

    print("=== 1-D Tumor kinetics, cycle residualize ON ===", flush=True)
    kin_on, _, _ = run_kinetics(adata, DEFAULT_GENES, residualize_v_cycle=True)
    print("=== 1-D Tumor kinetics, cycle residualize OFF ===", flush=True)
    kin_off, _, _ = run_kinetics(adata, DEFAULT_GENES, residualize_v_cycle=False)
    kin_on[["v_reox"]].join(kin_off["v_reox"], lsuffix="_cc", rsuffix="_nocc", how="inner").to_csv(
        OUT / "kinetics_v_reox_cycle_on_vs_off.csv"
    )
    r_cc, p_cc = spearman(kin_on["v_reox"].to_numpy(), kin_off["v_reox"].to_numpy())
    r_on_s, p_on_s = spearman(kin_on["v_reox"].to_numpy(), kin_on["cycle_s"].to_numpy())
    r_off_s, p_off_s = spearman(kin_off["v_reox"].to_numpy(), kin_off["cycle_s"].to_numpy())

    rows = [
        {
            "run": "tumor_1d_cc_on",
            "n": len(kin_on),
            "spearman_vHIF_vs_S": np.nan,
            "p_vHIF_vs_S": np.nan,
            "spearman_vCycle_vs_S": np.nan,
            "p_vCycle_vs_S": np.nan,
            "spearman_vHIF_vs_hifExpr": spearman(kin_on["v_reox"].to_numpy(), kin_on["hif_score"].to_numpy())[0],
            "p_vHIF_vs_hifExpr": spearman(kin_on["v_reox"].to_numpy(), kin_on["hif_score"].to_numpy())[1],
            "spearman_minus_vHIF_vs_our_v_reox": 1.0,
            "p_minus_vHIF_vs_our_v_reox": 0.0,
            "mean_abs_vCycle": np.nan,
            "mean_abs_vHIF": float(np.mean(np.abs(kin_on["v_reox"]))),
            "cycle_over_hif_ratio": np.nan,
            "frac_cycle_dominates": np.nan,
            "spearman_v_reox_vs_S": r_on_s,
            "p_v_reox_vs_S": p_on_s,
            "note": "our model; v residualized on E14 S/G2M",
        },
        {
            "run": "tumor_1d_cc_off",
            "n": len(kin_off),
            "spearman_vHIF_vs_hifExpr": spearman(kin_off["v_reox"].to_numpy(), kin_off["hif_score"].to_numpy())[0],
            "p_vHIF_vs_hifExpr": spearman(kin_off["v_reox"].to_numpy(), kin_off["hif_score"].to_numpy())[1],
            "spearman_v_reox_vs_S": r_off_s,
            "p_v_reox_vs_S": p_off_s,
            "spearman_v_reox_cc_on_vs_off": r_cc,
            "p_v_reox_cc_on_vs_off": p_cc,
            "mean_abs_vHIF": float(np.mean(np.abs(kin_off["v_reox"]))),
            "note": "our model; no v cycle residualization",
        },
    ]

    scv_runs = [
        ("scvelo_all_no_cc", "none"),
        ("scvelo_all_cc_allcells", "all"),
        ("scvelo_all_cc_e14fit", "e14"),
    ]
    per_cell = pd.DataFrame(
        {
            "sample": adata.obs["sample"].astype(str).values,
            "cell_group": adata.obs["cell_group"].astype(str).values,
            "cycle_s": s_score,
            "cycle_g2m": g2m_score,
            "hif_expr": adata.obs["hif_expr"].to_numpy(),
        },
        index=adata.obs_names,
    )
    per_cell = per_cell.join(kin_on[["v_reox", "theta_normoxic", "hypoxia_kinetics_state"]].rename(columns={"v_reox": "v_reox_cc_on"}))
    per_cell = per_cell.join(kin_off[["v_reox"]].rename(columns={"v_reox": "v_reox_cc_off"}))

    for name, mode in scv_runs:
        print(f"=== scVelo {name} ===", flush=True)
        b = run_scvelo(adata, mode, e14)
        v_hif = gene_set_velocity(b, HIF_THETA)
        v_cyc = gene_set_velocity(b, cycle_genes)
        aligned = pd.DataFrame({"v_hif": v_hif, "v_cyc": v_cyc}, index=b.obs_names)
        aligned = aligned.reindex(adata.obs_names)
        per_cell[f"{name}_vHIF"] = aligned["v_hif"].to_numpy()
        per_cell[f"{name}_vCycle"] = aligned["v_cyc"].to_numpy()
        ours = per_cell["v_reox_cc_on"].to_numpy()
        rows.append(
            summarize(
                name,
                aligned["v_hif"].to_numpy(),
                aligned["v_cyc"].to_numpy(),
                per_cell["hif_expr"].to_numpy(),
                s_score,
                ours,
            )
        )
        tumor = per_cell["cell_group"].eq("Tumor") & per_cell["v_reox_cc_on"].notna()
        r_t, p_t = spearman(-aligned.loc[tumor, "v_hif"].to_numpy(), per_cell.loc[tumor, "v_reox_cc_on"].to_numpy())
        rows[-1]["spearman_minus_vHIF_vs_our_v_reox_tumor"] = r_t
        rows[-1]["p_minus_vHIF_vs_our_v_reox_tumor"] = p_t
        rows[-1]["n_velocity_genes"] = int(b.var["velocity_genes"].sum()) if "velocity_genes" in b.var else np.nan
        print(f"  n_velocity_genes={rows[-1].get('n_velocity_genes')} tumor corr vs our v_reox r={r_t:.3f} p={p_t:.3g}", flush=True)

    summ = pd.DataFrame(rows)
    summ.to_csv(OUT / "scvelo_vs_kinetics_summary.csv", index=False)
    per_cell.to_csv(OUT / "scvelo_vs_kinetics_cells.csv")
    print("\n=== summary ===")
    print(summ.to_string(index=False))
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
