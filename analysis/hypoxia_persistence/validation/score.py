"""Score the E14/E15 HIF-down module on public human matrices."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.metrics import average_precision_score, roc_auc_score

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
GENES = ROOT / "reversion_key_genes.tsv"
QC = ROOT / "hypoxia_kinetics_gene_qc.csv"
THETA_MOUSE = [
    "Car9", "Vegfa", "Ldha", "Pdk1", "Egln1", "Pgk1", "Eno1", "Aldoa",
    "Angptl4", "P4ha1", "Ero1a", "Ankrd37", "Serpine1", "Cited2",
]


def human_theta_genes() -> pd.DataFrame:
    table = pd.read_csv(GENES, sep="\t")
    qc = pd.read_csv(QC)
    use = set(qc.loc[qc["use_theta"].eq(True), "gene"])
    table = table[table["mouse"].isin(use)].copy()
    w = qc.set_index("gene").loc[table["mouse"], "cohens_d_e15_vs_e14"].to_numpy(dtype=np.float64)
    w = np.clip(w, 0, None)
    w = w / w.sum()
    table = table.assign(weight=w)
    return table


def map_symbols(var_names: np.ndarray, symbols: list[str]) -> dict[str, int]:
    pos = {}
    for i, g in enumerate(var_names):
        s = str(g)
        if s.startswith("ENSG"):
            s = s.split(".")[0]
        pos.setdefault(s.upper(), i)
        if "_" in s:
            pos.setdefault(s.split("_")[0].upper(), i)
    return {s: pos[s.upper()] for s in symbols if s.upper() in pos}


def log1p_cp10k(X: np.ndarray) -> np.ndarray:
    lib = X.sum(axis=1)
    lib = np.where(lib <= 0, 1.0, lib)
    return np.log1p(X * (1e4 / lib)[:, None])


def module_score(X: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """X is cells x genes log1p; z-score genes then weighted mean (high = hypoxic)."""
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd = np.where(sd < 1e-12, 1.0, sd)
    z = (X - mu) / sd
    return z @ weights


def au_metrics(y_true: np.ndarray, score: np.ndarray) -> dict:
    y = y_true.astype(int)
    if y.min() == y.max():
        return {"auroc": np.nan, "auprc": np.nan, "n_pos": int(y.sum()), "n": int(y.size)}
    return {
        "auroc": float(roc_auc_score(y, score)),
        "auprc": float(average_precision_score(y, score)),
        "n_pos": int(y.sum()),
        "n": int(y.size),
    }


def summarize_split(
    score: np.ndarray,
    label: np.ndarray,
    pos_name: str,
    neg_name: str,
    *,
    robust: bool = False,
    n_boot: int = 400,
    n_perm: int = 400,
    seed: int = 0,
) -> dict:
    pos = label == pos_name
    neg = label == neg_name
    mask = pos | neg
    met = au_metrics(pos[mask], score[mask])
    d = np.nan
    if pos.sum() >= 5 and neg.sum() >= 5:
        v1, v2 = score[pos].var(ddof=1), score[neg].var(ddof=1)
        sp = np.sqrt(0.5 * (v1 + v2))
        d = float((score[pos].mean() - score[neg].mean()) / max(sp, 1e-12))
    met.update(
        {
            "pos": pos_name,
            "neg": neg_name,
            "median_pos": float(np.median(score[pos])) if pos.any() else np.nan,
            "median_neg": float(np.median(score[neg])) if neg.any() else np.nan,
            "cohens_d": d,
            "mwu_p": float(stats.mannwhitneyu(score[pos], score[neg], alternative="greater").pvalue)
            if pos.any() and neg.any()
            else np.nan,
        }
    )
    if robust and pos[mask].any() and neg[mask].any():
        y = pos[mask].astype(int)
        s = np.asarray(score)[mask]
        _, lo, hi = bootstrap_auroc(y, s, n=n_boot, seed=seed)
        _, p = perm_auroc(y, s, n=n_perm, seed=seed + 1)
        met["auroc_lo"] = lo
        met["auroc_hi"] = hi
        met["perm_p"] = p
    return met


def gene_auroc(X: np.ndarray, y: np.ndarray, symbols: list[str]) -> pd.DataFrame:
    rows = []
    for j, g in enumerate(symbols):
        met = au_metrics(y, X[:, j])
        rows.append({"gene": g, **met})
    return pd.DataFrame(rows)


def bootstrap_auroc(y: np.ndarray, score: np.ndarray, n: int = 400, seed: int = 0) -> tuple[float, float, float]:
    y = np.asarray(y).astype(int)
    score = np.asarray(score, dtype=np.float64)
    rng = np.random.default_rng(seed)
    vals = []
    nobs = y.size
    for _ in range(n):
        ix = rng.integers(0, nobs, nobs)
        yi, si = y[ix], score[ix]
        if yi.min() == yi.max():
            continue
        vals.append(float(roc_auc_score(yi, si)))
    if not vals:
        return np.nan, np.nan, np.nan
    arr = np.asarray(vals)
    return float(np.mean(arr)), float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))


def perm_auroc(y: np.ndarray, score: np.ndarray, n: int = 400, seed: int = 0) -> tuple[float, float]:
    y = np.asarray(y).astype(int)
    score = np.asarray(score, dtype=np.float64)
    obs = float(roc_auc_score(y, score)) if y.min() != y.max() else np.nan
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(n):
        yp = rng.permutation(y)
        if yp.min() == yp.max():
            continue
        if float(roc_auc_score(yp, score)) >= obs:
            ge += 1
    return obs, float((ge + 1) / (n + 1))


def perm_spearman(x: np.ndarray, y: np.ndarray, n: int = 400, seed: int = 0) -> tuple[float, float]:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    obs = float(stats.spearmanr(x, y).correlation)
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(n):
        if abs(float(stats.spearmanr(x, rng.permutation(y)).correlation)) >= abs(obs):
            ge += 1
    return obs, float((ge + 1) / (n + 1))


def visium_hex_moran(values, rows, cols) -> float:
    df = pd.DataFrame({"v": values, "r": rows, "c": cols}).dropna()
    if len(df) < 50:
        return np.nan
    key = {(int(r), int(c)): i for i, (r, c) in enumerate(zip(df.r, df.c))}
    v = df.v.to_numpy(dtype=np.float64)
    v = v - v.mean()
    num = 0.0
    wsum = 0.0
    neigh = [(-1, -1), (-1, 1), (0, -2), (0, 2), (1, -1), (1, 1)]
    for (r, c), i in key.items():
        for dr, dc in neigh:
            j = key.get((r + dr, c + dc))
            if j is not None:
                num += v[i] * v[j]
                wsum += 1
    den = float(np.sum(v * v))
    n = len(v)
    if wsum < 1 or den <= 0:
        return np.nan
    return float((n / wsum) * (num / den))


def perm_moran(values, rows, cols, n: int = 200, seed: int = 0) -> tuple[float, float]:
    obs = visium_hex_moran(values, rows, cols)
    if not np.isfinite(obs):
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    vals = np.asarray(values, dtype=np.float64)
    ge = 0
    for _ in range(n):
        if visium_hex_moran(rng.permutation(vals), rows, cols) >= obs:
            ge += 1
    return obs, float((ge + 1) / (n + 1))
