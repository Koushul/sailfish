#!/usr/bin/env python3
"""Library-level and HIF-gene spliced/unspliced QC for E14/E15."""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
DEFAULT_H5AD = Path("/ix1/ylee/kor11/MC38/E14SE15S/E14SE15S_gex_adt_placed.h5ad")
QC_GENES = HERE / "hypoxia_kinetics_gene_qc.csv"


def _summ(df: pd.DataFrame) -> pd.Series:
    return pd.Series(
        {
            "n": len(df),
            "med_spliced": float(df.spliced.median()),
            "med_unspliced": float(df.unspliced.median()),
            "med_ambiguous": float(df.ambiguous.median()),
            "med_u_over_s": float(df.u_over_s.median()),
            "med_frac_u": float(df.frac_u.median()),
            "med_frac_a": float(df.frac_a.median()),
            "med_genes_s": float(df.genes_s.median()),
            "med_genes_u": float(df.genes_u.median()),
            "p10_unspliced": float(df.unspliced.quantile(0.1)),
            "med_mapping": float(df.mapping_rate.median()),
        }
    )


def main() -> int:
    adata = ad.read_h5ad(DEFAULT_H5AD)
    S, U, Amb = adata.layers["spliced"].tocsr(), adata.layers["unspliced"].tocsr(), adata.layers["ambiguous"].tocsr()
    sample = adata.obs["sample"].astype(str).to_numpy()
    cg = adata.obs["cell_group"].astype(str).to_numpy()
    lib_s = np.asarray(S.sum(axis=1)).ravel()
    lib_u = np.asarray(U.sum(axis=1)).ravel()
    lib_a = np.asarray(Amb.sum(axis=1)).ravel()
    tot = lib_s + lib_u + lib_a
    cells = pd.DataFrame(
        {
            "sample": sample,
            "cell_group": cg,
            "spliced": lib_s,
            "unspliced": lib_u,
            "ambiguous": lib_a,
            "u_over_s": lib_u / np.maximum(lib_s, 1),
            "frac_u": lib_u / np.maximum(tot, 1),
            "frac_a": lib_a / np.maximum(tot, 1),
            "genes_s": np.asarray((S > 0).sum(axis=1)).ravel(),
            "genes_u": np.asarray((U > 0).sum(axis=1)).ravel(),
            "mapping_rate": adata.obs["mapping_rate"].to_numpy(),
            "tumor_qc": (cg == "Tumor") & (lib_s >= 5000),
        }
    )
    lib_rows = []
    for keys, g in cells.groupby(["sample", "cell_group"]):
        row = _summ(g)
        row["sample"], row["cell_group"], row["subset"] = keys[0], keys[1], "all"
        lib_rows.append(row)
    tqc = cells[cells.tumor_qc]
    for sm, g in tqc.groupby("sample"):
        row = _summ(g)
        row["sample"], row["cell_group"], row["subset"] = sm, "Tumor", "spliced_ge_5000"
        lib_rows.append(row)
    lib = pd.DataFrame(lib_rows)

    m15 = cells.tumor_qc & (cells["sample"] == "E15S")
    nz_u = np.asarray((U[np.flatnonzero(m15)] > 0).mean(axis=0)).ravel()
    mean_u = np.asarray(U[np.flatnonzero(m15)].mean(axis=0)).ravel()
    mean_s = np.asarray(S[np.flatnonzero(m15)].mean(axis=0)).ravel()
    lib.loc[lib.subset.eq("spliced_ge_5000") & lib["sample"].eq("E15S"), "n_genes_mean_u_ge1"] = int((mean_u >= 1).sum())
    lib.loc[lib.subset.eq("spliced_ge_5000") & lib["sample"].eq("E15S"), "n_genes_nz_u_ge_0.2"] = int((nz_u >= 0.2).sum())

    gqc = pd.read_csv(QC_GENES)
    pos = {str(g): i for i, g in enumerate(np.asarray(adata.var_names.astype(str)))}
    grow = []
    for _, r in gqc.iterrows():
        if r.gene not in pos:
            continue
        j = pos[r.gene]
        for sm in ["E14S", "E15S"]:
            m = cells.tumor_qc & (cells["sample"] == sm)
            idx = np.flatnonzero(m)
            s = np.asarray(S[idx, j].todense()).ravel()
            u = np.asarray(U[idx, j].todense()).ravel()
            am = np.asarray(Amb[idx, j].todense()).ravel()
            sp = float(stats.spearmanr(s, u).correlation) if s.std() > 0 and u.std() > 0 else np.nan
            grow.append(
                {
                    "gene": r.gene,
                    "use_theta": bool(r.use_theta),
                    "use_velocity": bool(r.use_velocity),
                    "sample": sm,
                    "n": int(m.sum()),
                    "mean_s": float(s.mean()),
                    "mean_u": float(u.mean()),
                    "mean_a": float(am.mean()),
                    "nz_s": float((s > 0).mean()),
                    "nz_u": float((u > 0).mean()),
                    "nz_a": float((am > 0).mean()),
                    "spearman_s_u": sp,
                    "kappa": float(u.mean() / (s.mean() + 1e-12)),
                }
            )
    genes = pd.DataFrame(grow)

    vgenes = gqc.loc[gqc.use_velocity.eq(True), "gene"]
    ix = [pos[g] for g in vgenes if g in pos]
    panel = []
    for sm in ["E14S", "E15S"]:
        m = cells.tumor_qc & (cells["sample"] == sm)
        idx = np.flatnonzero(m)
        s = np.asarray(S[idx][:, ix].sum(axis=1)).ravel()
        u = np.asarray(U[idx][:, ix].sum(axis=1)).ravel()
        panel.append(
            {
                "sample": sm,
                "n": int(m.sum()),
                "n_velocity_genes": len(ix),
                "median_panel_s": float(np.median(s)),
                "median_panel_u": float(np.median(u)),
                "p10_panel_u": float(np.quantile(u, 0.1)),
                "frac_panel_u_zero": float((u == 0).mean()),
            }
        )
    pan = pd.DataFrame(panel)

    lib.to_csv(HERE / "us_qc_library.csv", index=False)
    genes.to_csv(HERE / "us_qc_hif_genes.csv", index=False)
    pan.to_csv(HERE / "us_qc_velocity_panel.csv", index=False)
    print(lib.to_string(index=False))
    print(pan.to_string(index=False))
    print("wrote us_qc_library.csv us_qc_hif_genes.csv us_qc_velocity_panel.csv")
    print("E15 Tumor QC genes with mean unspliced >= 1:", int((mean_u >= 1).sum()), "of", U.shape[1])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
