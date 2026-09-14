#!/usr/bin/env python3
"""Persist vs revert: Tumor and neutrophils scored separately vs on one combined axis."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
NEU_PY = HERE.parent / "a223_neutrophil_hypoxia" / "run.py"
GENE_QC = HERE.parent / "a223_neutrophil_hypoxia" / "tumor_theta_gene_qc.csv"
E14E15 = Path("/ix1/ylee/kor11/MC38/E14SE15S/E14SE15S_gex_adt_placed.h5ad")
OUT = Path("/ix1/ylee/kor11/A223/tumor_kinetics/lineage_combine")
REPO = HERE / "results"
TUMOR_UMI = 5000.0
NEU_UMI = 1500.0
T_LOW, T_HIGH = 0.30, 0.70


def load_neu():
    spec = importlib.util.spec_from_file_location("a223_neu_hif", NEU_PY)
    neu = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(neu)
    return neu


def logistic_theta(h: np.ndarray, h_lo: float, h_hi: float) -> np.ndarray:
    span = max(h_hi - h_lo, 1e-6)
    h0 = 0.5 * (h_lo + h_hi)
    tau = span / 6.0
    return 1.0 / (1.0 + np.exp((h - h0) / tau))


def state(theta: np.ndarray) -> np.ndarray:
    out = np.full(theta.shape, "partial", dtype=object)
    out[theta <= T_LOW] = "persistent"
    out[theta >= T_HIGH] = "reverted"
    return out


def calibrate(log_s: np.ndarray, w: np.ndarray, ref: np.ndarray, hyp: np.ndarray) -> dict:
    mu = log_s[ref].mean(axis=0)
    sd = log_s[ref].std(axis=0)
    sd = np.where(sd < 1e-12, 1.0, sd)
    h = ((log_s - mu) / sd) @ w
    h_lo = float(np.median(h[ref]))
    q = float(np.quantile(h[hyp], 0.75))
    hi = hyp & (h >= q)
    h_hi = float(np.median(h[hi])) if hi.any() else float(np.median(h[hyp]))
    return {"mu": mu, "sd": sd, "h_lo": h_lo, "h_hi": h_hi, "n_ref": int(ref.sum()), "n_hyp": int(hyp.sum())}


def apply_cal(log_s: np.ndarray, w: np.ndarray, cal: dict) -> np.ndarray:
    h = ((log_s - cal["mu"]) / cal["sd"]) @ w
    return logistic_theta(h, cal["h_lo"], cal["h_hi"])


def row(cohort: str, lineage: str, calibration: str, theta: np.ndarray) -> dict:
    st = state(theta)
    n = int(theta.size)
    return {
        "cohort": cohort,
        "lineage": lineage,
        "calibration": calibration,
        "n": n,
        "persistent_pct": 100.0 * (st == "persistent").mean() if n else np.nan,
        "partial_pct": 100.0 * (st == "partial").mean() if n else np.nan,
        "reverted_pct": 100.0 * (st == "reverted").mean() if n else np.nan,
        "median_theta": float(np.median(theta)) if n else np.nan,
    }


def plot_compare(df: pd.DataFrame, cohort: str, path: Path) -> None:
    sub = df[df["cohort"] == cohort].copy()
    sub["label"] = sub["lineage"] + "\n" + sub["calibration"]
    fig, ax = plt.subplots(figsize=(9.2, 4.6))
    x = np.arange(len(sub))
    persist = sub["persistent_pct"].to_numpy()
    partial = sub["partial_pct"].to_numpy()
    reverted = sub["reverted_pct"].to_numpy()
    ax.bar(x, persist, color="#c0392b", label="persistent")
    ax.bar(x, partial, bottom=persist, color="#f4d03f", label="partial")
    ax.bar(x, reverted, bottom=persist + partial, color="#2980b9", label="reverted")
    ax.set_xticks(x)
    ax.set_xticklabels(sub["label"].tolist(), rotation=25, ha="right", fontsize=8)
    ax.set_ylim(0, 100)
    ax.set_ylabel("% QC cells")
    ax.set_title(f"{cohort}: separate vs combined neutrophil axis")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    REPO.mkdir(parents=True, exist_ok=True)
    neu = load_neu()
    mapping = neu.load_symbols()
    qc = pd.read_csv(GENE_QC)
    use = qc["use_theta"].astype(str).str.lower().isin(["true", "1"])
    genes = qc.loc[use, "gene"].astype(str).tolist()
    w = qc.set_index("gene").loc[genes, "cohens_d_e15_vs_e14"].clip(lower=0).to_numpy(dtype=np.float64)
    w = w / w.sum()

    ref = ad.read_h5ad(E14E15)
    neu.set_symbols(ref, mapping)
    lib_r = np.asarray(ref.layers["spliced"].sum(axis=1), dtype=np.float64).ravel()
    sample = ref.obs["sample"].astype(str).to_numpy()
    group = ref.obs["cell_group"].astype(str).to_numpy()
    e14_t = (group == "Tumor") & (lib_r >= TUMOR_UMI) & (sample == "E14S")
    e15_t = (group == "Tumor") & (lib_r >= TUMOR_UMI) & (sample == "E15S")
    e14_n = (group == "neutrophil") & (lib_r >= NEU_UMI) & (sample == "E14S")
    e15_n = (group == "neutrophil") & (lib_r >= NEU_UMI) & (sample == "E15S")
    e14_j = e14_t | e14_n
    e15_j = e15_t | e15_n
    log_r = np.log1p(
        neu.size_normalize(
            neu.csr_cols(ref.layers["spliced"], [ {str(g): i for i, g in enumerate(ref.var_names.astype(str))}[g] for g in genes ]),
            lib_r,
        )
    )
    cal_t = calibrate(log_r, w, e14_t, e15_t)
    cal_n = calibrate(log_r, w, e14_n, e15_n)
    cal_j = calibrate(log_r, w, e14_j, e15_j)
    th_r_t = apply_cal(log_r, w, cal_t)
    th_r_n = apply_cal(log_r, w, cal_n)
    th_r_j = apply_cal(log_r, w, cal_j)

    a223 = neu.load_merged(mapping)
    scores = neu.annotate(a223)
    lib_a = scores["spliced_umi"].to_numpy()
    pos = {str(g): i for i, g in enumerate(a223.var_names.astype(str))}
    log_a = np.log1p(neu.size_normalize(neu.csr_cols(a223.layers["spliced"], [pos[g] for g in genes]), lib_a))
    tmask = scores["palak_tumor"].to_numpy() & (lib_a >= TUMOR_UMI)
    nmask = scores["palak_neutrophil"].to_numpy() & (lib_a >= NEU_UMI)
    both = tmask | nmask
    th_a_t = apply_cal(log_a, w, cal_t)
    th_a_n = apply_cal(log_a, w, cal_n)
    th_a_j = apply_cal(log_a, w, cal_j)

    rows = [
        row("E15S_MC38", "Tumor", "separate (Tumor E14)", th_r_t[e15_t]),
        row("E15S_MC38", "neutrophil", "separate (neu E14)", th_r_n[e15_n]),
        row("E15S_MC38", "Tumor", "combined (Tumor+neu E14)", th_r_j[e15_t]),
        row("E15S_MC38", "neutrophil", "combined (Tumor+neu E14)", th_r_j[e15_n]),
        row("E15S_MC38", "Tumor+neutrophil", "pooled on Tumor axis", np.concatenate([th_r_t[e15_t], th_r_t[e15_n]])),
        row("E15S_MC38", "Tumor+neutrophil", "pooled on neu axis", np.concatenate([th_r_n[e15_t], th_r_n[e15_n]])),
        row("E15S_MC38", "Tumor+neutrophil", "pooled on combined axis", th_r_j[e15_j]),
        row("E15S_MC38", "Tumor+neutrophil", "separate then concat states", np.concatenate([th_r_t[e15_t], th_r_n[e15_n]])),
        row("A223_E27E29", "Tumor", "separate (Tumor E14)", th_a_t[tmask]),
        row("A223_E27E29", "neutrophil", "separate (neu E14)", th_a_n[nmask]),
        row("A223_E27E29", "Tumor", "combined (Tumor+neu E14)", th_a_j[tmask]),
        row("A223_E27E29", "neutrophil", "combined (Tumor+neu E14)", th_a_j[nmask]),
        row("A223_E27E29", "Tumor+neutrophil", "pooled on Tumor axis", np.concatenate([th_a_t[tmask], th_a_t[nmask]])),
        row("A223_E27E29", "Tumor+neutrophil", "pooled on neu axis", np.concatenate([th_a_n[tmask], th_a_n[nmask]])),
        row("A223_E27E29", "Tumor+neutrophil", "pooled on combined axis", th_a_j[both]),
        row("A223_E27E29", "Tumor+neutrophil", "separate then concat states", np.concatenate([th_a_t[tmask], th_a_n[nmask]])),
    ]
    table = pd.DataFrame(rows)
    cal_meta = pd.DataFrame(
        [
            {"calibration": "Tumor", **{k: cal_t[k] for k in ("h_lo", "h_hi", "n_ref", "n_hyp")}},
            {"calibration": "neutrophil", **{k: cal_n[k] for k in ("h_lo", "h_hi", "n_ref", "n_hyp")}},
            {"calibration": "combined", **{k: cal_j[k] for k in ("h_lo", "h_hi", "n_ref", "n_hyp")}},
        ]
    )

    table.to_csv(OUT / "persist_separate_vs_combined.csv", index=False)
    table.to_csv(REPO / "persist_separate_vs_combined.csv", index=False)
    cal_meta.to_csv(OUT / "calibration_anchors.csv", index=False)
    cal_meta.to_csv(REPO / "lineage_calibration_anchors.csv", index=False)
    plot_compare(table, "E15S_MC38", OUT / "e15_separate_vs_combined.png")
    plot_compare(table, "A223_E27E29", OUT / "a223_separate_vs_combined.png")
    print(cal_meta.to_string(index=False))
    print(table.to_string(index=False))
    del ref, a223


if __name__ == "__main__":
    main()
