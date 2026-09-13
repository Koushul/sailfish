#!/usr/bin/env python3
"""Test whether an Image-IT / GFP-channel sort would explain TAM in E15S.

Image-IT LIVE Green ROS (carboxy-H2DCFDA) reports in the GFP/FITC band.
TAMs autofluoresce there (flavins, lipofuscin) and also make real ROS (NOX2).
This script compares E14 vs E15 composition and TAM autofluor vs ROS programs.
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
DEFAULT_H5AD = Path("/ix1/ylee/kor11/MC38/E14SE15S/E14SE15S_gex_adt_placed.h5ad")

SETS = {
    "tam": ("C1qa", "C1qb", "C1qc", "Adgre1", "Csf1r", "Cd68", "Fcgr1", "Trem2", "Mrc1", "Mertk"),
    "lyso": ("Lamp1", "Ctsb", "Ctsd", "Hexb", "Cd63", "Ctsl", "Psap"),
    "nox": ("Cybb", "Ncf1", "Ncf2", "Ncf4", "Nox1", "Nox4"),
    "nrf2": ("Hmox1", "Nqo1", "Sod2", "Txnrd1", "Gclm", "Gclc"),
    "hif": ("Ldha", "Vegfa", "Pgk1", "P4ha1", "Bnip3", "Slc2a1"),
    "mono": ("Ly6c2", "Ccr2", "Chil3"),
}


def log1p_cp10k_sum(S, lib, ix):
    if not ix:
        return np.zeros(S.shape[0])
    x = np.asarray(S[:, ix].sum(axis=1)).ravel()
    return np.log1p(x * (1e4 / np.maximum(lib, 1)))


def fisher_or_e15(sample, is_g):
    y15 = sample == "E15S"
    tab = np.array(
        [
            [(~y15 & ~is_g).sum(), (~y15 & is_g).sum()],
            [(y15 & ~is_g).sum(), (y15 & is_g).sum()],
        ]
    )
    or_, p = stats.fisher_exact(tab)
    return int(tab[0, 1]), int(tab[0].sum()), int(tab[1, 1]), int(tab[1].sum()), float(or_), float(p)


def main() -> int:
    a = ad.read_h5ad(DEFAULT_H5AD)
    var = np.asarray(a.var_names.astype(str))
    S = a.layers["spliced"].tocsr()
    lib = np.asarray(S.sum(axis=1)).ravel()
    mt = np.array([str(g).startswith("mt-") for g in var])
    frac_mt = np.asarray(S[:, mt].sum(axis=1)).ravel() / np.maximum(lib, 1)
    pos = {g: i for i, g in enumerate(var)}
    sample = a.obs["sample"].astype(str).to_numpy()
    cg = a.obs["cell_group"].astype(str).to_numpy()
    df = pd.DataFrame({"sample": sample, "cell_group": cg, "lib": lib, "frac_mt": frac_mt})
    for name, genes in SETS.items():
        df[name] = log1p_cp10k_sum(S, lib, [pos[g] for g in genes if g in pos])

    comp = []
    for g in ["TAM", "Tumor", "neutrophil", "others"]:
        n14, N14, n15, N15, or_, p = fisher_or_e15(sample, cg == g)
        comp.append(
            {
                "cell_group": g,
                "n_E14": n14,
                "share_E14": n14 / N14,
                "n_E15": n15,
                "share_E15": n15 / N15,
                "OR_E15_vs_E14": or_,
                "fisher_p": p,
            }
        )
    cdf = pd.DataFrame(comp)

    feat_rows = []
    for group in ["TAM", "neutrophil", "Tumor"]:
        sub = df[df.cell_group == group]
        for col in ["lib", "frac_mt", "tam", "lyso", "nox", "nrf2", "hif", "mono"]:
            if col not in sub:
                continue
            a14 = sub.loc[sub["sample"] == "E14S", col].to_numpy()
            a15 = sub.loc[sub["sample"] == "E15S", col].to_numpy()
            if min(len(a14), len(a15)) < 20:
                continue
            d = (a15.mean() - a14.mean()) / np.sqrt(0.5 * (a14.var(ddof=1) + a15.var(ddof=1)))
            feat_rows.append(
                {
                    "cell_group": group,
                    "feature": col,
                    "median_E14": float(np.median(a14)),
                    "median_E15": float(np.median(a15)),
                    "cohens_d_E15_minus_E14": float(d),
                    "mwu_p": float(stats.mannwhitneyu(a15, a14).pvalue),
                    "auroc_E15": float(roc_auc_score(np.r_[np.zeros(len(a14)), np.ones(len(a15))], np.r_[a14, a15])),
                }
            )
    fdf = pd.DataFrame(feat_rows)
    cdf.to_csv(HERE / "ros_sort_composition.csv", index=False)
    fdf.to_csv(HERE / "ros_sort_programs.csv", index=False)
    print(cdf.to_string(index=False))
    print(fdf.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
