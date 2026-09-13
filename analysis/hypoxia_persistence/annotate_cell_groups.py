#!/usr/bin/env python3
"""Four-way cell groups on E14SE15S_gex_adt_placed.h5ad: Tumor, TAM, neutrophil, others.

Protein ADT is too sparse to use. Groups are assigned from RNA clusters
(PCA of HVGs + k-means) using Ptprc, C1q/Adgre1, and S100a8/Cxcr2.
Monocytes, T/NK, DC, and low-quality cells are others.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

DEFAULT_H5AD = Path("/ix1/ylee/kor11/MC38/E14SE15S/E14SE15S_gex_adt_placed.h5ad")
HERE = Path(__file__).resolve().parent

TAM_GENES = ("C1qa", "C1qb", "C1qc", "Adgre1", "Csf1r", "Cd68", "Fcgr1", "Trem2", "Mrc1", "Mertk")
NEU_GENES = ("S100a8", "S100a9", "Cxcr2", "Csf3r", "Retnlg", "Lcn2")
MONO_GENES = ("Ly6c2", "Ccr2", "Chil3")
TNK_GENES = ("Cd3d", "Cd3e", "Cd8a", "Nkg7")
DC_GENES = ("Flt3", "Xcr1", "Clec9a")


def log_cp10k(X):
    X = X.astype(np.float64)
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    lib = np.asarray(X.sum(axis=1)).ravel()
    lib[lib <= 0] = 1.0
    X = X.multiply(1e4 / lib[:, None])
    X.data = np.log1p(X.data)
    return X.tocsr()


def hvg_index(X, n_hvg: int = 2000, min_cells: int = 30) -> np.ndarray:
    nnz = np.asarray((X > 0).sum(axis=0)).ravel()
    idx = np.where(nnz >= min_cells)[0]
    mean = np.asarray(X[:, idx].mean(axis=0)).ravel()
    mean_sq = np.asarray(X[:, idx].multiply(X[:, idx]).mean(axis=0)).ravel()
    var = mean_sq - mean**2
    n_hvg = min(n_hvg, idx.size)
    return idx[np.argsort(var)[-n_hvg:]]


def gene_mean(X, var_names: np.ndarray, genes: tuple[str, ...]) -> np.ndarray:
    want = set(genes)
    cols = [i for i, g in enumerate(var_names) if g in want]
    if not cols:
        return np.zeros(X.shape[0], dtype=np.float64)
    return np.asarray(X[:, cols].mean(axis=1)).ravel()


def cluster_map(med: pd.DataFrame) -> dict[int, str]:
    """Label each cluster from median marker scores."""
    out = {}
    for k, r in med.iterrows():
        if r["ptprc"] < 0.4:
            out[int(k)] = "others" if r["gex"] < 2000 else "Tumor"
            continue
        if r["neu"] >= 1.2 and r["neu"] > r["tam"] + 0.3:
            out[int(k)] = "neutrophil"
            continue
        if r["tam"] >= 0.9 and r["tam"] >= r["mono"] - 0.15 and r["tnk"] < 0.5:
            out[int(k)] = "TAM"
            continue
        out[int(k)] = "others"
    return out


def annotate(adata, n_clusters: int = 18, n_pcs: int = 30, seed: int = 0) -> pd.DataFrame:
    var_names = np.asarray(adata.var_names.astype(str))
    X = log_cp10k(adata.layers["spliced"] if "spliced" in adata.layers else adata.X)
    hvg = hvg_index(X)
    z = StandardScaler().fit_transform(X[:, hvg].toarray())
    pcs = PCA(n_components=n_pcs, random_state=seed).fit_transform(z)
    cl = MiniBatchKMeans(n_clusters=n_clusters, random_state=seed, batch_size=2048, n_init=10).fit_predict(pcs)

    scores = pd.DataFrame(
        {
            "cluster": cl,
            "ptprc": gene_mean(X, var_names, ("Ptprc",)),
            "tam": gene_mean(X, var_names, TAM_GENES),
            "neu": gene_mean(X, var_names, NEU_GENES),
            "mono": gene_mean(X, var_names, MONO_GENES),
            "tnk": gene_mean(X, var_names, TNK_GENES),
            "dc": gene_mean(X, var_names, DC_GENES),
            "gex": np.asarray(adata.obs["gex_counts"], dtype=np.float64),
        },
        index=adata.obs_names,
    )
    med = scores.groupby("cluster")[["ptprc", "tam", "neu", "mono", "tnk", "dc", "gex"]].median()
    mapping = cluster_map(med)
    scores["cell_group"] = scores["cluster"].map(mapping)

    neu_cell = (scores["neu"] > scores["tam"] + 0.4) & (scores["neu"] > 1.5) & (scores["ptprc"] > 0.4)
    tam_cell = (scores["tam"] > 1.2) & (scores["tam"] > scores["neu"] + 0.3) & (scores["tnk"] < 0.6) & (scores["ptprc"] > 0.4)
    tumor_cell = (scores["ptprc"] < 0.25) & (scores["gex"] >= 2000) & (scores["tam"] < 0.8) & (scores["neu"] < 1.0)
    scores.loc[neu_cell, "cell_group"] = "neutrophil"
    scores.loc[tam_cell & ~neu_cell, "cell_group"] = "TAM"
    scores.loc[tumor_cell & (scores["cell_group"] == "others"), "cell_group"] = "Tumor"

    scores["sample"] = adata.obs["sample"].astype(str).values
    if "barcodes" in adata.obs:
        scores["barcodes"] = adata.obs["barcodes"].astype(str).values
    print("cluster medians:\n", med.round(3).to_string())
    print("cluster → group:", {k: mapping[k] for k in sorted(mapping)})
    print(scores.groupby(["sample", "cell_group"]).size().unstack(fill_value=0))
    print(scores["cell_group"].value_counts().to_string())
    return scores


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", default=str(DEFAULT_H5AD))
    ap.add_argument("--out-csv", default=str(HERE / "cell_groups.csv"))
    ap.add_argument("--write-h5ad", action="store_true", help="write cell_group into the h5ad obs")
    args = ap.parse_args()

    import anndata as ad

    path = Path(args.h5ad)
    adata = ad.read_h5ad(path)
    scores = annotate(adata)
    out = Path(args.out_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    scores.to_csv(out)
    print(f"wrote {out}")
    if args.write_h5ad:
        adata.obs["cell_group"] = scores["cell_group"].astype("category")
        adata.obs["cell_group_cluster"] = scores["cluster"].astype("int32")
        adata.write_h5ad(path)
        print(f"updated obs on {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
