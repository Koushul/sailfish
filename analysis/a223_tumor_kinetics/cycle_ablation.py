#!/usr/bin/env python3
"""Tumor 1-D θ with vs without cell cycle (A223 + E15)."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
NEU = HERE.parent / "a223_neutrophil_hypoxia" / "run.py"
GENE_QC = HERE.parent / "a223_neutrophil_hypoxia" / "tumor_theta_gene_qc.csv"
CYCLE_JSON = HERE / "tirosh_cell_cycle_mouse.json"
E14E15 = Path("/ix1/ylee/kor11/MC38/E14SE15S/E14SE15S_gex_adt_placed.h5ad")
OUT = Path("/ix1/ylee/kor11/A223/tumor_kinetics/cycle_ablation")
REPO = HERE / "results"
TUMOR_UMI = 5000.0
T_LOW, T_HIGH = 0.30, 0.70
R2_DROP = 0.10


def load_neu():
    spec = importlib.util.spec_from_file_location("a223_neu_hif", NEU)
    neu = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(neu)
    return neu


def mean_log1p(norm: np.ndarray, names: np.ndarray, symbols: list[str]) -> np.ndarray:
    pos = {str(g): i for i, g in enumerate(names)}
    cols = [pos[g] for g in symbols if g in pos]
    if len(cols) < 5:
        raise SystemExit(f"cycle score: {len(cols)} genes")
    return np.log1p(norm[:, cols]).mean(axis=1)


def ols_cycle(y: np.ndarray, s: np.ndarray, g2m: np.ndarray, mask: np.ndarray, r2_min: float = 0.05):
    y14 = y[mask]
    X = np.column_stack([np.ones(int(mask.sum())), s[mask], g2m[mask]])
    coef, *_ = np.linalg.lstsq(X, y14, rcond=None)
    pred = X @ coef
    sst = np.sum((y14 - y14.mean()) ** 2)
    r2 = 1.0 - np.sum((y14 - pred) ** 2) / sst if sst > 0 else 0.0
    if r2 < r2_min:
        return 0.0, 0.0, float(r2)
    return float(coef[1]), float(coef[2]), float(r2)


def logistic_theta(h: np.ndarray, h_lo: float, h_hi: float) -> np.ndarray:
    span = max(h_hi - h_lo, 1e-6)
    h0 = 0.5 * (h_lo + h_hi)
    tau = span / 6.0
    return 1.0 / (1.0 + np.exp((h - h0) / tau))


def theta_from_log(log_s: np.ndarray, w: np.ndarray, e14: np.ndarray, e15: np.ndarray) -> np.ndarray:
    mu = log_s[e14].mean(axis=0)
    sd = log_s[e14].std(axis=0)
    sd = np.where(sd < 1e-12, 1.0, sd)
    h = ((log_s - mu) / sd) @ w
    h_lo = float(np.median(h[e14]))
    q = float(np.quantile(h[e15], 0.75))
    hi = e15 & (h >= q)
    h_hi = float(np.median(h[hi])) if hi.any() else float(np.median(h[e15]))
    return logistic_theta(h, h_lo, h_hi)


def state(theta: np.ndarray) -> np.ndarray:
    out = np.full(theta.shape, "partial", dtype=object)
    out[theta <= T_LOW] = "persistent"
    out[theta >= T_HIGH] = "reverted"
    return out


def frac_row(name: str, theta: np.ndarray, s: np.ndarray, g2m: np.ndarray) -> dict:
    st = state(theta)
    n = theta.size
    rho_s = spearmanr(theta, s)
    rho_g = spearmanr(theta, g2m)
    return {
        "model": name,
        "n": int(n),
        "persistent_pct": 100.0 * (st == "persistent").mean(),
        "partial_pct": 100.0 * (st == "partial").mean(),
        "reverted_pct": 100.0 * (st == "reverted").mean(),
        "median_theta": float(np.median(theta)),
        "spearman_theta_vs_S": float(rho_s.statistic),
        "p_theta_vs_S": float(rho_s.pvalue),
        "spearman_theta_vs_G2M": float(rho_g.statistic),
        "p_theta_vs_G2M": float(rho_g.pvalue),
    }


def residualize(log_s: np.ndarray, s: np.ndarray, g2m: np.ndarray, fit_mask: np.ndarray, s0: float, g0: float):
    adj = np.array(log_s, copy=True)
    rows = []
    for j in range(log_s.shape[1]):
        bs, bg, r2 = ols_cycle(log_s[:, j], s, g2m, fit_mask)
        adj[:, j] = log_s[:, j] - bs * (s - s0) - bg * (g2m - g0)
        rows.append({"j": j, "b_S": bs, "b_G2M": bg, "cycle_r2": r2})
    return adj, pd.DataFrame(rows)


def cycle_scores(neu, adata, lib, cycle) -> tuple[np.ndarray, np.ndarray]:
    var_names = np.asarray(adata.var_names.astype(str))
    pos = {str(g): i for i, g in enumerate(var_names)}
    syms = [g for g in cycle["S"] + cycle["G2M"] if g in pos]
    cols = [pos[g] for g in syms]
    names = np.array(syms)
    norm = neu.size_normalize(neu.csr_cols(adata.layers["spliced"], cols), lib)
    s = mean_log1p(norm, names, [g for g in cycle["S"] if g in pos])
    g2m = mean_log1p(norm, names, [g for g in cycle["G2M"] if g in pos])
    return s, g2m


def hif_log(neu, adata, lib, genes: list[str]) -> np.ndarray:
    pos = {str(g): i for i, g in enumerate(adata.var_names.astype(str))}
    cols = [pos[g] for g in genes]
    return np.log1p(neu.size_normalize(neu.csr_cols(adata.layers["spliced"], cols), lib))


def plot_models(rows: pd.DataFrame, title: str, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    x = np.arange(len(rows))
    persist = rows["persistent_pct"].to_numpy()
    partial = rows["partial_pct"].to_numpy()
    reverted = rows["reverted_pct"].to_numpy()
    ax.bar(x, persist, color="#c0392b", label="persistent")
    ax.bar(x, partial, bottom=persist, color="#f4d03f", label="partial")
    ax.bar(x, reverted, bottom=persist + partial, color="#2980b9", label="reverted")
    ax.set_xticks(x)
    ax.set_xticklabels(rows["model"].tolist(), rotation=20, ha="right")
    ax.set_ylim(0, 100)
    ax.set_ylabel("% QC Tumor")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_theta_vs_s(theta_off, theta_on, s, labels, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.8), sharey=True)
    for ax, th, lab in zip(axes, (theta_off, theta_on), labels):
        ax.scatter(s, th, s=4, alpha=0.15, c="#34495e", linewidths=0)
        ax.set_xlabel("Tirosh S score")
        ax.set_title(lab)
        ax.axhline(T_LOW, color="#c0392b", ls="--", lw=0.8)
        ax.axhline(T_HIGH, color="#2980b9", ls="--", lw=0.8)
    axes[0].set_ylabel("θ_normoxic")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    REPO.mkdir(parents=True, exist_ok=True)
    neu = load_neu()
    mapping = neu.load_symbols()
    cycle = json.loads(CYCLE_JSON.read_text())
    qc = pd.read_csv(GENE_QC)
    use = qc["use_theta"].astype(str).str.lower().isin(["true", "1"])
    genes_all = qc.loc[use, "gene"].astype(str).tolist()
    w_all = qc.set_index("gene").loc[genes_all, "cohens_d_e15_vs_e14"].clip(lower=0).to_numpy(dtype=np.float64)
    w_all = w_all / w_all.sum()
    keep_drop = qc.set_index("gene").loc[genes_all, "cycle_r2_spliced_e14"].lt(R2_DROP).to_numpy()
    genes_drop = [g for g, k in zip(genes_all, keep_drop) if k]
    w_drop = qc.set_index("gene").loc[genes_drop, "cohens_d_e15_vs_e14"].clip(lower=0).to_numpy(dtype=np.float64)
    w_drop = w_drop / w_drop.sum()

    ref = ad.read_h5ad(E14E15)
    neu.set_symbols(ref, mapping)
    lib_r = np.asarray(ref.layers["spliced"].sum(axis=1), dtype=np.float64).ravel()
    sample = ref.obs["sample"].astype(str).to_numpy()
    group = ref.obs["cell_group"].astype(str).to_numpy()
    tumor_r = (group == "Tumor") & (lib_r >= TUMOR_UMI)
    e14 = tumor_r & (sample == "E14S")
    e15 = tumor_r & (sample == "E15S")
    s_r, g2m_r = cycle_scores(neu, ref, lib_r, cycle)
    log_r = hif_log(neu, ref, lib_r, genes_all)
    log_r_drop = hif_log(neu, ref, lib_r, genes_drop)
    s0, g0 = float(s_r[e14].mean()), float(g2m_r[e14].mean())
    log_r_cc, gene_cycle = residualize(log_r, s_r, g2m_r, e14, s0, g0)
    gene_cycle.insert(0, "gene", genes_all)

    theta_e15_off = theta_from_log(log_r, w_all, e14, e15)[e15]
    theta_e15_cc = theta_from_log(log_r_cc, w_all, e14, e15)[e15]
    theta_e15_drop = theta_from_log(log_r_drop, w_drop, e14, e15)[e15]
    s_e15, g_e15 = s_r[e15], g2m_r[e15]

    e15_tbl = pd.DataFrame(
        [
            frac_row("no_cycle (default θ)", theta_e15_off, s_e15, g_e15),
            frac_row("residualize S+G2M on HIF genes", theta_e15_cc, s_e15, g_e15),
            frac_row(f"drop genes with cycle R²≥{R2_DROP}", theta_e15_drop, s_e15, g_e15),
        ]
    )
    e15_tbl.insert(0, "cohort", "E15S_MC38")

    a223 = neu.load_merged(mapping)
    scores = neu.annotate(a223)
    lib_a = scores["spliced_umi"].to_numpy()
    palak_t = scores["palak_tumor"].to_numpy() & (lib_a >= TUMOR_UMI)
    s_a, g2m_a = cycle_scores(neu, a223, lib_a, cycle)
    log_a = hif_log(neu, a223, lib_a, genes_all)
    log_a_drop = hif_log(neu, a223, lib_a, genes_drop)
    log_a_cc, _ = residualize(log_a, s_a, g2m_a, palak_t, float(s_a[palak_t].mean()), float(g2m_a[palak_t].mean()))

    def score_a223(log_ref, log_tgt, w):
        mu = log_ref[e14].mean(axis=0)
        sd = log_ref[e14].std(axis=0)
        sd = np.where(sd < 1e-12, 1.0, sd)
        h_ref = ((log_ref - mu) / sd) @ w
        h_lo = float(np.median(h_ref[e14]))
        q = float(np.quantile(h_ref[e15], 0.75))
        hi = e15 & (h_ref >= q)
        h_hi = float(np.median(h_ref[hi])) if hi.any() else float(np.median(h_ref[e15]))
        h = ((log_tgt - mu) / sd) @ w
        return logistic_theta(h, h_lo, h_hi)

    theta_a_off = score_a223(log_r, log_a, w_all)[palak_t]
    theta_a_cc = score_a223(log_r_cc, log_a_cc, w_all)[palak_t]
    theta_a_drop = score_a223(log_r_drop, log_a_drop, w_drop)[palak_t]
    s_at, g_at = s_a[palak_t], g2m_a[palak_t]
    subtype = scores.loc[palak_t, "palak_cell_type"].astype(str).to_numpy()

    a_tbl = pd.DataFrame(
        [
            frac_row("no_cycle (default θ)", theta_a_off, s_at, g_at),
            frac_row("residualize S+G2M on HIF genes", theta_a_cc, s_at, g_at),
            frac_row(f"drop genes with cycle R²≥{R2_DROP}", theta_a_drop, s_at, g_at),
        ]
    )
    a_tbl.insert(0, "cohort", "A223_E27E29")

    compare = pd.concat([e15_tbl, a_tbl], ignore_index=True)

    sub_rows = []
    for model, th in (
        ("no_cycle", theta_a_off),
        ("residualize", theta_a_cc),
        ("drop_cycle_genes", theta_a_drop),
    ):
        df = pd.DataFrame({"theta": th, "subtype": subtype, "S": s_at})
        df["state"] = state(df["theta"].to_numpy())
        for sub, g in df.groupby("subtype"):
            sub_rows.append(
                {
                    "model": model,
                    "palak_cell_type": sub,
                    "n": int(len(g)),
                    "persistent_pct": 100.0 * (g["state"] == "persistent").mean(),
                    "reverted_pct": 100.0 * (g["state"] == "reverted").mean(),
                    "median_theta": float(g["theta"].median()),
                    "median_S": float(g["S"].median()),
                }
            )
    subtype_tbl = pd.DataFrame(sub_rows)

    agree = pd.crosstab(
        pd.Series(state(theta_a_off), name="no_cycle"),
        pd.Series(state(theta_a_cc), name="residualize"),
        normalize=False,
    )
    agree_drop = pd.crosstab(
        pd.Series(state(theta_a_off), name="no_cycle"),
        pd.Series(state(theta_a_drop), name="drop_cycle_genes"),
        normalize=False,
    )

    gene_cycle.to_csv(OUT / "hif_gene_cycle_ols_e14.csv", index=False)
    compare.to_csv(OUT / "persist_by_cycle_model.csv", index=False)
    subtype_tbl.to_csv(OUT / "persist_by_subtype_and_model.csv", index=False)
    agree.to_csv(OUT / "a223_crosstab_no_cycle_vs_residualize.csv")
    agree_drop.to_csv(OUT / "a223_crosstab_no_cycle_vs_drop.csv")
    compare.to_csv(REPO / "persist_by_cycle_model.csv", index=False)
    subtype_tbl.to_csv(REPO / "persist_by_subtype_and_model.csv", index=False)
    gene_cycle.to_csv(REPO / "hif_gene_cycle_ols_e14.csv", index=False)

    plot_models(e15_tbl, "E15S Tumor θ with vs without cell cycle", OUT / "e15_stacked_cycle.png")
    plot_models(a_tbl, "A223 Tumor θ with vs without cell cycle", OUT / "a223_stacked_cycle.png")
    plot_theta_vs_s(
        theta_a_off,
        theta_a_cc,
        s_at,
        ["A223 no cycle correction", "A223 residualize S+G2M"],
        OUT / "a223_theta_vs_S.png",
    )
    dropped = [g for g, k in zip(genes_all, keep_drop) if not k]
    print("dropped genes", dropped)
    print(compare.to_string(index=False))
    print(subtype_tbl.to_string(index=False))
    print("label agreement residualize\n", agree)
    print("label agreement drop\n", agree_drop)
    del ref, a223


if __name__ == "__main__":
    main()
