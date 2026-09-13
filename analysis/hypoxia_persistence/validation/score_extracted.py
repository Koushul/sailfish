#!/usr/bin/env python3
"""Score extracted Seurat dumps + write a combined validation summary."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse, stats

from score import gene_auroc, human_theta_genes, log1p_cp10k, module_score, summarize_split

HERE = Path(__file__).resolve().parent
EXT = HERE / "data" / "extracted"
OUT = HERE / "results"

HASH_GSE200207 = {
    "Hashtag1": "normal_Nx",
    "Hashtag2": "normal_Hx",
    "Hashtag3": "tumor_Nx",
    "Hashtag4": "tumor_Nx",
    "Hashtag5": "normal_Nx_lowserum",
}

TIME_GSE296547 = {
    "NOR": 0,
    "HYPO6": 6,
    "HYPO15": 15,
    "HYPO30": 30,
    "HYPO30NOR6": 36,
    "HYPO30NOR15": 45,
    "HYPO30NOR30": 60,
}


def load_extracted(prefix: Path):
    genes = pd.read_csv(str(prefix) + "_genes.tsv", header=None)[0].astype(str).tolist()
    cells = pd.read_csv(str(prefix) + "_cells.tsv", header=None)[0].astype(str).tolist()
    meta = pd.read_csv(str(prefix) + "_meta.csv", index_col=0)
    M = mmread(str(prefix) + "_counts.mtx")
    if sparse.issparse(M):
        M = M.tocsr()
    X = np.asarray(M.T.todense() if M.shape[0] == len(genes) else M.todense(), dtype=np.float64)
    if X.shape[0] != len(cells):
        X = X.T
    meta = meta.reindex(cells)
    return np.array(genes), np.array(cells), X, meta


def aligned_module(gene_names, X):
    table = human_theta_genes()
    ix, w = [], []
    used = []
    up = {g.upper(): i for i, g in enumerate(gene_names)}
    for g, wt in zip(table["human"], table["weight"]):
        if g.upper() in up:
            ix.append(up[g.upper()])
            w.append(wt)
            used.append(g)
    w = np.asarray(w, dtype=np.float64)
    w = w / w.sum()
    logx = log1p_cp10k(X[:, ix])
    return module_score(logx, w), used, logx, w


def moran_i(values, rows, cols):
    from score import visium_hex_moran

    return visium_hex_moran(values, rows, cols)


def gse200207():
    genes, cells, X, meta = load_extracted(EXT / "GSE200207")
    score, used, logx, w = aligned_module(genes, X)
    lab = meta["hash.ID"].map(HASH_GSE200207).fillna("other").to_numpy()
    rows = []
    met = summarize_split(score, lab, "normal_Hx", "normal_Nx")
    met.update({"dataset": "GSE200207", "split": "normal kidney 0.5% vs 21% O2 (module)", "genes": ",".join(used)})
    rows.append(met)
    y = (lab == "normal_Hx")
    mask = np.isin(lab, ["normal_Hx", "normal_Nx"])
    gtab = gene_auroc(logx[mask], y[mask], used)
    gtab.to_csv(OUT / "GSE200207_gene_auroc_normal.csv", index=False)
    rows.append(
        {
            "dataset": "GSE200207",
            "split": "normal module vs median gene",
            "auroc": met["auroc"],
            "median_gene_auroc": float(gtab.auroc.median()),
            "best_gene": gtab.sort_values("auroc").iloc[-1]["gene"],
            "best_gene_auroc": float(gtab.auroc.max()),
            "n": met["n"],
        }
    )
    met2 = summarize_split(score, np.where(np.isin(lab, ["tumor_Nx"]), "tumor_Nx", lab), "tumor_Nx", "normal_Nx")
    met2.update({"dataset": "GSE200207", "split": "ccRCC tumor 21% vs normal 21% (constitutive HIF)"})
    rows.append(met2)
    pd.DataFrame({"cell": cells, "label": lab, "hif_module": score}).to_csv(OUT / "GSE200207_cell_scores.csv", index=False)
    return pd.DataFrame(rows)


def gse296547():
    genes, cells, X, meta = load_extracted(EXT / "GSE296547")
    score, used, logx, w = aligned_module(genes, X)
    tp = meta["timepoint"].astype(str).to_numpy()
    rows = []
    # hypoxia vs never-hypoxic
    lab = np.where(np.isin(tp, ["HYPO6", "HYPO15", "HYPO30"]), "Hx", np.where(tp == "NOR", "Nx", "reox"))
    met = summarize_split(score, lab, "Hx", "Nx")
    met.update({"dataset": "GSE296547", "split": "HYPO 6-30d vs NOR (module)"})
    rows.append(met)
    # reox vs still hypoxic (HYPO30)
    lab2 = np.where(np.isin(tp, ["HYPO30NOR6", "HYPO30NOR15", "HYPO30NOR30"]), "reox", np.where(tp == "HYPO30", "Hx", "other"))
    met2 = summarize_split(-score, lab2, "reox", "Hx")  # reox should have LOWER hypoxia score
    # rewrite with hypoxia score: pos=Hx vs reox
    met2 = summarize_split(score, lab2, "Hx", "reox")
    met2.update({"dataset": "GSE296547", "split": "HYPO30 vs reoxygenation (module; Hx should stay higher)"})
    rows.append(met2)
    # monotonicity: hypoxia days
    order = ["NOR", "HYPO6", "HYPO15", "HYPO30", "HYPO30NOR6", "HYPO30NOR15", "HYPO30NOR30"]
    med = {t: float(np.median(score[tp == t])) for t in order if (tp == t).any()}
    rho = float(stats.spearmanr([TIME_GSE296547[t] for t in order if t in med and not t.startswith("HYPO30NOR")], [med[t] for t in order if t in med and not t.startswith("HYPO30NOR")]).correlation)
    rows.append({"dataset": "GSE296547", "split": "spearman module vs hypoxia day (NOR+HYPO only)", "spearman_hours": rho, **{f"median_{k}": v for k, v in med.items()}})
    pd.DataFrame({"cell": cells, "timepoint": tp, "hif_module": score, "cell_type": meta["cell_type"].astype(str).values}).to_csv(
        OUT / "GSE296547_cell_scores.csv", index=False
    )
    return pd.DataFrame(rows)


def visium_moran():
    rows = []
    for prefix in ["GSM7688224_Tumor-897", "GSM7688225_Tumor-899"]:
        p = OUT / f"{prefix}_spot_scores.csv"
        if not p.exists():
            continue
        df = pd.read_csv(p)
        in_t = df["in_tissue"].eq(1) if "in_tissue" in df.columns else np.ones(len(df), bool)
        mi = moran_i(df.loc[in_t, "hif_module"].to_numpy(), df.loc[in_t, "row"], df.loc[in_t, "col"])
        rows.append({"dataset": "GSE240212", "split": f"{prefix} Moran I of HIF module", "moran_i": mi, "n": int(in_t.sum())})
    return pd.DataFrame(rows)


def main():
    from robustness import main as robust_main

    return robust_main()


if __name__ == "__main__":
    main()
