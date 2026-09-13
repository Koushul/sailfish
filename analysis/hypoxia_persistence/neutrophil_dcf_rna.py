#!/usr/bin/env python3
"""RNA proxies for DCF+/primed neutrophils (burst itself is post-transcriptional)."""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
DEFAULT_H5AD = Path("/ix1/ylee/kor11/MC38/E14SE15S/E14SE15S_gex_adt_placed.h5ad")
NEU = HERE / "neutrophil_kinetics.csv"

# Capacity (not burst): NADPH oxidase subunits. Burst is p47 phosphorylation / assembly.
NOX = ("Cybb", "Cyba", "Ncf1", "Ncf2", "Ncf4", "Rac2")
# Wright 2013 PLoS ONE: 1 h TNF/GM-CSF priming common transcripts (not CYBB).
PRIMING = ("Cxcl1", "Cxcl2", "Il1a", "Il1b", "Il1rn", "Icam1")
PRIMING_TNF = ("Ccl3", "Ccl4", "Tnf", "Nfkbia", "Tnfaip3")
NRF2 = ("Hmox1", "Nqo1", "Gclc", "Gclm", "Txnrd1", "Sod2", "Gpx1")
HIF = ("Ldha", "Vegfa", "Pgk1", "Eno1", "Slc2a1", "P4ha1")
MATURE = ("Fpr1", "Cxcr2", "C5ar1")


def module(S, lib, pos, genes):
    ix = [pos[g] for g in genes if g in pos]
    if not ix:
        return np.zeros(S.shape[0]), []
    x = np.asarray(S[:, ix].sum(axis=1)).ravel()
    return np.log1p(x * (1e4 / np.maximum(lib, 1))), [g for g in genes if g in pos]


def contrast(score, pos_mask, neg_mask):
    a, b = score[pos_mask], score[neg_mask]
    if min(a.size, b.size) < 10:
        return {}
    d = (a.mean() - b.mean()) / np.sqrt(0.5 * (a.var(ddof=1) + b.var(ddof=1)))
    y = np.r_[np.ones(a.size), np.zeros(b.size)]
    s = np.r_[a, b]
    return {
        "n_pos": int(a.size),
        "n_neg": int(b.size),
        "median_pos": float(np.median(a)),
        "median_neg": float(np.median(b)),
        "cohens_d": float(d),
        "mwu_p": float(stats.mannwhitneyu(a, b, alternative="greater").pvalue),
        "auroc": float(roc_auc_score(y, s)),
    }


def main() -> int:
    a = ad.read_h5ad(DEFAULT_H5AD)
    var = np.asarray(a.var_names.astype(str))
    pos = {g: i for i, g in enumerate(var)}
    S = a.layers["spliced"].tocsr()
    lib = np.asarray(S.sum(axis=1)).ravel()
    sample = a.obs["sample"].astype(str).to_numpy()
    cg = a.obs["cell_group"].astype(str).to_numpy()
    neu = cg == "neutrophil"
    kin = pd.read_csv(NEU, index_col=0)
    state = pd.Series(np.array(["unlabeled"] * a.n_obs, dtype=object), index=a.obs_names)
    state.loc[kin.index] = kin["hypoxia_kinetics_state"].astype(str).values
    state = state.to_numpy()

    sets = {
        "nox_capacity": NOX,
        "priming_wright2013": PRIMING,
        "priming_tnf": PRIMING_TNF,
        "nrf2_ros_response": NRF2,
        "hif_tumor_module": HIF,
        "mature_receptors": MATURE,
    }
    rows = []
    gene_rows = []
    for name, genes in sets.items():
        sc, used = module(S, lib, pos, genes)
        for g in used:
            j = pos[g]
            raw = np.asarray(S[:, j].todense()).ravel()
            m14 = neu & (sample == "E14S")
            m15 = neu & (sample == "E15S")
            gene_rows.append(
                {
                    "module": name,
                    "gene": g,
                    "nz_neu_E14": float((raw[m14] > 0).mean()),
                    "nz_neu_E15": float((raw[m15] > 0).mean()),
                    "mean_log1p_E14": float(np.log1p(raw[m14] * (1e4 / np.maximum(lib[m14], 1))).mean()),
                    "mean_log1p_E15": float(np.log1p(raw[m15] * (1e4 / np.maximum(lib[m15], 1))).mean()),
                }
            )
        met = contrast(sc, neu & (sample == "E15S"), neu & (sample == "E14S"))
        met.update({"split": "E15 vs E14 neutrophils", "module": name, "genes": ",".join(used)})
        rows.append(met)
        e15p = neu & (sample == "E15S") & (state == "persistent")
        e15r = neu & (sample == "E15S") & (state == "reverted")
        met2 = contrast(sc, e15p, e15r)
        met2.update({"split": "E15 persistent vs reverted neutrophils", "module": name, "genes": ",".join(used)})
        rows.append(met2)
    pd.DataFrame(rows).to_csv(HERE / "neutrophil_dcf_rna.csv", index=False)
    pd.DataFrame(gene_rows).to_csv(HERE / "neutrophil_dcf_rna_genes.csv", index=False)
    print(pd.DataFrame(rows)[["split", "module", "auroc", "cohens_d", "mwu_p", "median_pos", "median_neg"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
