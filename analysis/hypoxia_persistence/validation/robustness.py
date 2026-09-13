#!/usr/bin/env python3
"""Recompute labeled splits with bootstrap CIs, label permutation, LOO, and Visium Moran I."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from score import (
    gene_auroc,
    human_theta_genes,
    perm_moran,
    perm_spearman,
    summarize_split,
)
from score_extracted import TIME_GSE296547, aligned_module, load_extracted

HERE = Path(__file__).resolve().parent
EXT = HERE / "data" / "extracted"
OUT = HERE / "results"


def _loo_module(logx: np.ndarray, weights: np.ndarray, used: list[str], y: np.ndarray) -> pd.DataFrame:
    from score import module_score

    rows = []
    w = np.asarray(weights, dtype=np.float64)
    for j, g in enumerate(used):
        keep = np.ones(len(used), dtype=bool)
        keep[j] = False
        ww = w[keep]
        ww = ww / ww.sum()
        score = module_score(logx[:, keep], ww)
        met = gene_auroc(score[:, None], y, [f"drop_{g}"]).iloc[0].to_dict()
        met["dropped"] = g
        rows.append(met)
    return pd.DataFrame(rows)


def gse200207() -> tuple[pd.DataFrame, pd.DataFrame]:
    from score_extracted import HASH_GSE200207

    genes, cells, X, meta = load_extracted(EXT / "GSE200207")
    score, used, logx, w = aligned_module(genes, X)
    lab = meta["hash.ID"].map(HASH_GSE200207).fillna("other").to_numpy()
    rows = []
    met = summarize_split(score, lab, "normal_Hx", "normal_Nx", robust=True, seed=1)
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
            "worst_gene": gtab.sort_values("auroc").iloc[0]["gene"],
            "worst_gene_auroc": float(gtab.auroc.min()),
            "n": met["n"],
        }
    )
    loo = _loo_module(logx[mask], w, used, y[mask].astype(int))
    loo.to_csv(OUT / "GSE200207_loo_normal.csv", index=False)
    rows.append(
        {
            "dataset": "GSE200207",
            "split": "normal LOO min AUROC (drop one gene)",
            "auroc": float(loo.auroc.min()),
            "worst_gene": loo.sort_values("auroc").iloc[0]["dropped"],
            "n": met["n"],
        }
    )
    met2 = summarize_split(
        score,
        np.where(np.isin(lab, ["tumor_Nx"]), "tumor_Nx", lab),
        "tumor_Nx",
        "normal_Nx",
        robust=True,
        seed=2,
    )
    met2.update({"dataset": "GSE200207", "split": "ccRCC tumor 21% vs normal 21% (constitutive HIF)"})
    rows.append(met2)
    pd.DataFrame({"cell": cells, "label": lab, "hif_module": score}).to_csv(OUT / "GSE200207_cell_scores.csv", index=False)
    return pd.DataFrame(rows), loo


def gse296547() -> pd.DataFrame:
    genes, cells, X, meta = load_extracted(EXT / "GSE296547")
    score, used, logx, w = aligned_module(genes, X)
    tp = meta["timepoint"].astype(str).to_numpy()
    ct = meta["cell_type"].astype(str).to_numpy()
    at2 = np.isin(ct, ["Hypoxic_AT2_cells", "Normoxic_AT2_cells"])
    rows = []
    lab = np.where(np.isin(tp, ["HYPO6", "HYPO15", "HYPO30"]), "Hx", np.where(tp == "NOR", "Nx", "reox"))
    met = summarize_split(score, lab, "Hx", "Nx", robust=True, seed=3)
    met.update({"dataset": "GSE296547", "split": "all cells HYPO 6-30d vs NOR"})
    rows.append(met)
    met_at2 = summarize_split(score[at2], lab[at2], "Hx", "Nx", robust=True, seed=4)
    met_at2.update({"dataset": "GSE296547", "split": "AT2 only HYPO 6-30d vs NOR (avoids hypoxic-cluster circularity)"})
    rows.append(met_at2)
    lab2 = np.where(np.isin(tp, ["HYPO30NOR6", "HYPO30NOR15", "HYPO30NOR30"]), "reox", np.where(tp == "HYPO30", "Hx", "other"))
    met2 = summarize_split(score, lab2, "Hx", "reox", robust=True, seed=5)
    met2.update({"dataset": "GSE296547", "split": "all cells HYPO30 vs reoxygenation"})
    rows.append(met2)
    met2a = summarize_split(score[at2], lab2[at2], "Hx", "reox", robust=True, seed=6)
    met2a.update({"dataset": "GSE296547", "split": "AT2 only HYPO30 vs reoxygenation"})
    rows.append(met2a)
    order_hx = ["NOR", "HYPO6", "HYPO15", "HYPO30"]
    med = {t: float(np.median(score[tp == t])) for t in TIME_GSE296547 if (tp == t).any()}
    med_at2 = {t: float(np.median(score[at2 & (tp == t)])) for t in TIME_GSE296547 if (at2 & (tp == t)).any()}
    rho, p = perm_spearman(
        np.array([TIME_GSE296547[t] for t in order_hx]),
        np.array([med[t] for t in order_hx]),
        n=400,
        seed=7,
    )
    rows.append(
        {
            "dataset": "GSE296547",
            "split": "spearman module vs hypoxia day (all cells, NOR+HYPO)",
            "spearman_hours": rho,
            "perm_p": p,
            **{f"median_{k}": v for k, v in med.items()},
        }
    )
    rho2, p2 = perm_spearman(
        np.array([TIME_GSE296547[t] for t in order_hx if t in med_at2]),
        np.array([med_at2[t] for t in order_hx if t in med_at2]),
        n=400,
        seed=8,
    )
    rows.append(
        {
            "dataset": "GSE296547",
            "split": "spearman module vs hypoxia day (AT2, NOR+HYPO)",
            "spearman_hours": rho2,
            "perm_p": p2,
            **{f"median_{k}": v for k, v in med_at2.items()},
        }
    )
    pd.DataFrame({"cell": cells, "timepoint": tp, "hif_module": score, "cell_type": ct, "at2": at2}).to_csv(
        OUT / "GSE296547_cell_scores.csv", index=False
    )
    return pd.DataFrame(rows)


def gse292771() -> pd.DataFrame:
    df = pd.read_csv(OUT / "GSE292771_cell_scores.csv")
    genes = [g for g in human_theta_genes()["human"] if g in df.columns]
    rows = []
    for i, geno in enumerate(["WT", "HIF1KO", "HIF2KO"]):
        sub = df[df.genotype == geno]
        met = summarize_split(sub.hif_module.to_numpy(), sub.oxygen.to_numpy(), "Hx", "Nx", robust=True, seed=10 + i)
        met.update({"dataset": "GSE292771", "split": f"{geno} Hx vs Nx (module)"})
        rows.append(met)
        y = (sub.oxygen == "Hx").to_numpy()
        gtab = gene_auroc(sub[genes].to_numpy(), y, genes)
        gtab.to_csv(OUT / f"GSE292771_gene_auroc_{geno}.csv", index=False)
    wt = pd.read_csv(OUT / "GSE292771_gene_auroc_WT.csv")
    rows.append(
        {
            "dataset": "GSE292771",
            "split": "WT module vs median gene AUROC",
            "auroc": rows[0]["auroc"],
            "median_gene_auroc": float(wt.auroc.median()),
            "best_gene": wt.sort_values("auroc").iloc[-1]["gene"],
            "best_gene_auroc": float(wt.auroc.max()),
            "n": rows[0]["n"],
            "n_pos": rows[0]["n_pos"],
        }
    )
    wt_df = df[df.genotype == "WT"]
    X = wt_df[genes].to_numpy()
    y = (wt_df.oxygen == "Hx").to_numpy().astype(int)
    w = human_theta_genes().set_index("human").loc[genes, "weight"].to_numpy()
    loo = _loo_module(X, w, genes, y)
    loo.to_csv(OUT / "GSE292771_loo_WT.csv", index=False)
    rows.append(
        {
            "dataset": "GSE292771",
            "split": "WT LOO min AUROC (drop one gene)",
            "auroc": float(loo.auroc.min()),
            "worst_gene": loo.sort_values("auroc").iloc[0]["dropped"],
            "n": rows[0]["n"],
        }
    )
    return pd.DataFrame(rows)


def gse30019() -> pd.DataFrame:
    per = pd.read_csv(OUT / "GSE30019_sample_scores.csv")
    score = per.hif_module.to_numpy()
    hours = per.hours.to_numpy()
    lab = np.where(hours == 0, "Hx", np.where(hours == 24, "Nx", "mid"))
    met = summarize_split(score, lab, "Hx", "Nx", robust=True, n_boot=200, n_perm=200, seed=20)
    rho, p = perm_spearman(hours, score, n=400, seed=21)
    met.update({"dataset": "GSE30019", "split": "0h hypoxia vs 24h reox (bulk module)", "spearman_hours": rho, "perm_p_spearman": p})
    genes = pd.read_csv(OUT / "GSE30019_gene_reox.csv")
    met["best_drop_gene"] = genes.sort_values("log2fc_24_vs_0").iloc[0]["gene"]
    met["best_drop_log2fc"] = float(genes["log2fc_24_vs_0"].min())
    return pd.DataFrame([met])


def gse240212() -> pd.DataFrame:
    rows = []
    for i, prefix in enumerate(["GSM7688224_Tumor-897", "GSM7688225_Tumor-899"]):
        df = pd.read_csv(OUT / f"{prefix}_spot_scores.csv")
        in_t = df["in_tissue"].eq(1) if "in_tissue" in df.columns else np.ones(len(df), bool)
        s = df.loc[in_t, "hif_module"].to_numpy()
        r = df.loc[in_t, "row"].to_numpy()
        c = df.loc[in_t, "col"].to_numpy()
        mi, p = perm_moran(s, r, c, n=200, seed=30 + i)
        rows.append(
            {
                "dataset": "GSE240212",
                "split": f"{prefix} Visium hex Moran I of HIF module",
                "moran_i": mi,
                "perm_p": p,
                "n": int(in_t.sum()),
                "median_module": float(np.median(s)),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    OUT.mkdir(exist_ok=True)
    parts = []
    print("=== GSE200207 ===")
    t, _ = gse200207()
    print(t[["dataset", "split", "auroc", "auroc_lo", "auroc_hi", "perm_p"]].to_string(index=False))
    parts.append(t)
    print("=== GSE296547 ===")
    t = gse296547()
    print(t[["dataset", "split", "auroc", "auroc_lo", "auroc_hi", "perm_p", "spearman_hours"]].to_string(index=False))
    parts.append(t)
    print("=== GSE292771 ===")
    t = gse292771()
    print(t[["dataset", "split", "auroc", "auroc_lo", "auroc_hi", "perm_p"]].to_string(index=False))
    parts.append(t)
    print("=== GSE30019 ===")
    t = gse30019()
    print(t.to_string(index=False))
    parts.append(t)
    print("=== GSE240212 ===")
    t = gse240212()
    print(t.to_string(index=False))
    parts.append(t)
    out = pd.concat(parts, ignore_index=True, sort=False)
    out.to_csv(OUT / "summary.csv", index=False)
    print("wrote", OUT / "summary.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
