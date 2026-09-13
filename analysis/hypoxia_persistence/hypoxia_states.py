#!/usr/bin/env python3
"""Hypoxia persistence / reversion from E14 (normoxic) vs E15 (hypoxic) scRNA-seq.

GFP+ in HIF fate-mapping marks *history* of hypoxia (protein persists after
reoxygenation). The transcriptome can still look hypoxic (persistent) or
normoxic (reverted). This module scores HIF/hypoxia programs, then classifies
E15 cells relative to the E14 baseline.

Without a GFP column, every E15 cell is treated as hypoxia-exposed; pass
--gfp-obs to restrict reversion calls to reporter-positive cells.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

HERE = Path(__file__).resolve().parent
GENESET_DIR = HERE / "genesets"
DEFAULT_H5AD = Path("/ix1/ylee/kor11/MC38/E14SE15S/E14SE15S_gex_adt_placed.h5ad")
DEFAULT_ID2NAME = Path("/ix1/ylee/Palak/MC38/genesets/gene_id_to_name.tsv")
DEFAULT_WAGNER = Path("/ix1/ylee/shared/external/data/WagnerCollab/mc38_velocity.h5ad")
NORMOXIA_SAMPLE = "E14S"
HYPOXIA_SAMPLE = "E15S"


def load_id2name(path: Path) -> dict[str, str]:
    df = pd.read_csv(path, sep="\t", header=None, names=["gene_id", "symbol"])
    return dict(zip(df["gene_id"].astype(str), df["symbol"].astype(str)))


def geneset_symbols(name: str) -> list[str]:
    if name == "hif1a_mouse":
        df = pd.read_csv(GENESET_DIR / "hif1a_mouse_ensembl.tsv", sep="\t")
        return df["gene_name"].astype(str).tolist()
    p = GENESET_DIR / f"{name}.json"
    payload = json.loads(p.read_text())
    key = next(iter(payload))
    return list(payload[key]["geneSymbols"])


def resolve_var_index(var_names: np.ndarray, id2name: dict[str, str], symbols: list[str]) -> np.ndarray:
    """Match gene-set symbols to AnnData var (Ensembl ids or symbols), case-insensitive."""
    want = {s.lower() for s in symbols}
    name_l = {gid: sym.lower() for gid, sym in id2name.items()}
    idx = []
    for i, v in enumerate(var_names):
        vs = str(v)
        if vs.lower() in want:
            idx.append(i)
            continue
        if name_l.get(vs, "") in want:
            idx.append(i)
    return np.asarray(idx, dtype=int)


def log1p_mean(X, cols: np.ndarray) -> np.ndarray:
    sub = X[:, cols]
    if sparse.issparse(sub):
        sub = sub.astype(np.float64)
        sub.data = np.log1p(sub.data)
        return np.asarray(sub.mean(axis=1)).ravel()
    return np.log1p(np.asarray(sub, dtype=np.float64)).mean(axis=1)


def score_genesets(adata, id2name: dict[str, str]) -> pd.DataFrame:
    sets = {
        "hallmark_hypoxia": "HALLMARK_HYPOXIA",
        "buffa_hypoxia": "BUFFA_HYPOXIA_METAGENE",
        "winter_hypoxia": "WINTER_HYPOXIA_METAGENE",
        "hif1a_targets": "hif1a_mouse",
    }
    vn = np.asarray(adata.var_names.astype(str))
    out = {}
    n_hit = {}
    for col, gs in sets.items():
        idx = resolve_var_index(vn, id2name, geneset_symbols(gs))
        n_hit[col] = int(idx.size)
        if idx.size < 5:
            raise SystemExit(f"{col}: only {idx.size} genes matched")
        out[col] = log1p_mean(adata.X, idx)
    print("geneset hits:", n_hit)
    return pd.DataFrame(out, index=adata.obs_names)


def classify(
    obs: pd.DataFrame,
    scores: pd.DataFrame,
    gfp: np.ndarray | None,
    gfp_min: float,
    n_gmm: int,
) -> pd.DataFrame:
    df = obs[["sample"]].copy()
    df = df.join(scores)
    axis = df["hallmark_hypoxia"].to_numpy(dtype=np.float64)
    e14 = (df["sample"].astype(str) == NORMOXIA_SAMPLE).to_numpy()
    e15 = (df["sample"].astype(str) == HYPOXIA_SAMPLE).to_numpy()
    if e14.sum() < 20 or e15.sum() < 20:
        raise SystemExit(f"need both samples; got E14={e14.sum()} E15={e15.sum()}")

    # Extreme-contrast classifier: E14 vs E15 cells in the top quartile of the hypoxia axis.
    e15_ax = axis[e15]
    hi = np.quantile(e15_ax, 0.75)
    y = np.zeros(df.shape[0], dtype=int)
    pos = e15 & (axis >= hi)
    y[pos] = 1
    train = e14 | pos
    X = StandardScaler().fit_transform(scores.to_numpy(dtype=np.float64))
    clf = LogisticRegression(max_iter=500)
    clf.fit(X[train], y[train])
    p_hyp = clf.predict_proba(X)[:, 1]
    df["P_hypoxic"] = p_hyp
    df["P_normoxic"] = 1.0 - p_hyp

    gmm = GaussianMixture(n_components=n_gmm, covariance_type="full", random_state=0)
    gmm.fit(axis[e15].reshape(-1, 1))
    order = np.argsort(gmm.means_.ravel())  # low → high hypoxia
    raw = gmm.predict(axis.reshape(-1, 1))
    remap = {int(old): int(new) for new, old in enumerate(order)}
    df["gmm_component"] = np.array([remap[int(k)] for k in raw], dtype=int)
    df["gmm_P_high"] = gmm.predict_proba(axis.reshape(-1, 1))[:, int(order[-1])]

    gfp_pos = np.ones(df.shape[0], dtype=bool) if gfp is None else (np.asarray(gfp) >= gfp_min)
    df["gfp_positive"] = gfp_pos
    df["gfp_available"] = gfp is not None

    state = np.array(["unassigned"] * df.shape[0], dtype=object)
    state[e14] = "never_hypoxic"
    # E15: GMM high = persistent transcriptome; GMM low = reverted; middle = partial
    e15_comp = df["gmm_component"].to_numpy()
    if n_gmm == 2:
        state[e15 & (e15_comp == 1)] = "persistent"
        state[e15 & (e15_comp == 0)] = "reverted"
    else:
        state[e15 & (e15_comp == n_gmm - 1)] = "persistent"
        state[e15 & (e15_comp == 0)] = "reverted"
        state[e15 & (e15_comp > 0) & (e15_comp < n_gmm - 1)] = "partial"
    if gfp is not None:
        # Reporter history: GFP− E15 cells are not called reverted (never labeled).
        state[e15 & ~gfp_pos & (state == "reverted")] = "unlabeled_normoxic"
        state[e15 & ~gfp_pos & (state == "persistent")] = "unlabeled_hypoxic"
    df["hypoxia_state"] = state
    df["persistence_score"] = np.where(e15, p_hyp, np.nan)
    df["reversion_score"] = np.where(e15, 1.0 - p_hyp, np.nan)
    if gfp is not None:
        z_gfp = (np.asarray(gfp, dtype=np.float64) - np.nanmean(gfp)) / (np.nanstd(gfp) + 1e-9)
        z_rna = (axis - axis.mean()) / (axis.std() + 1e-9)
        df["gfp_rna_discordance"] = z_gfp - z_rna
    return df


def compare_wagner(df: pd.DataFrame, wagner_h5ad: Path) -> pd.DataFrame | None:
    if not wagner_h5ad.is_file():
        return None
    import anndata as ad

    w = ad.read_h5ad(wagner_h5ad, backed="r")
    wobs = w.obs[["barcodes", "sample", "hypoxia_state", "hypoxia_state_gmm", "P_normoxic", "hypoxia_score"]].copy()
    wobs["join"] = wobs["barcodes"].astype(str) + "-" + wobs["sample"].astype(str)
    left = df.reset_index().rename(columns={"index": "obs_name"})
    if "barcodes" in left.columns:
        left["join"] = left["barcodes"].astype(str) + "-" + left["sample"].astype(str)
    else:
        left["join"] = left["obs_name"].astype(str)
    m = left.merge(wobs, on="join", how="inner", suffixes=("", "_w"))
    print(f"Wagner barcode overlap: {len(m)} / {len(df)}")
    if m.empty:
        return None
    ct = pd.crosstab(m["hypoxia_state"], m["hypoxia_state_w"])
    print("new state vs Wagner hypoxia_state:\n", ct)
    return m


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", default=str(DEFAULT_H5AD))
    ap.add_argument("--id2name", default=str(DEFAULT_ID2NAME))
    ap.add_argument("--gfp-obs", default=None, help="obs column with GFP (or HIF reporter) intensity")
    ap.add_argument("--gfp-min", type=float, default=None, help="threshold for GFP+ (default: median of nonzero)")
    ap.add_argument("--n-gmm", type=int, default=3, choices=(2, 3))
    ap.add_argument("--wagner", default=str(DEFAULT_WAGNER), help="optional prior labels for overlap table")
    ap.add_argument("--out", default=str(HERE / "hypoxia_states.csv"))
    args = ap.parse_args()

    import anndata as ad

    adata = ad.read_h5ad(args.h5ad)
    id2name = load_id2name(Path(args.id2name))
    scores = score_genesets(adata, id2name)
    gfp = None
    gfp_min = 0.0
    if args.gfp_obs:
        if args.gfp_obs not in adata.obs.columns:
            raise SystemExit(f"--gfp-obs {args.gfp_obs} not in obs")
        gfp = adata.obs[args.gfp_obs].to_numpy(dtype=np.float64)
        gfp_min = args.gfp_min if args.gfp_min is not None else float(np.median(gfp[gfp > 0])) if np.any(gfp > 0) else 0.0
        print(f"GFP+ threshold {gfp_min:.4g}  n+={(gfp >= gfp_min).sum()}")
    df = classify(adata.obs, scores, gfp, gfp_min, args.n_gmm)
    if "barcodes" in adata.obs:
        df["barcodes"] = adata.obs["barcodes"].astype(str).values
    for c in ("map_row", "map_col"):
        if c in adata.obs:
            df[c] = adata.obs[c].values
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out)
    print(f"wrote {out} n={len(df)}")
    print(df.groupby(["sample", "hypoxia_state"]).size().unstack(fill_value=0))
    print(df.groupby("hypoxia_state")[["hallmark_hypoxia", "P_hypoxic"]].median())
    compare_wagner(df, Path(args.wagner))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
