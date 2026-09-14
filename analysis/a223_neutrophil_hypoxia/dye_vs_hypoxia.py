#!/usr/bin/env python3
"""Test: DCF+/FITC+ neutrophils from dye use vs true hypoxia."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
NEU_PY = HERE / "run.py"
CELLS = Path("/ix1/ylee/kor11/A223/neutrophil_hypoxia/cells.csv")
E14E15 = Path("/ix1/ylee/kor11/MC38/E14SE15S/E14SE15S_gex_adt_placed.h5ad")
GENE_QC = HERE / "tumor_theta_gene_qc.csv"
OUT = Path("/ix1/ylee/kor11/A223/neutrophil_hypoxia/dye_vs_hypoxia")
REPO = HERE
NEU_UMI = 1500.0
TUMOR_UMI = 5000.0
T_LOW, T_HIGH = 0.30, 0.70
MPO_EST = ("Mpo", "Ces1d", "Ces1g", "Ces2a", "Ces2c", "Ces2e")
HIF = ("Ldha", "Vegfa", "Pgk1", "Eno1", "Slc2a1", "P4ha1")


def load_neu():
    spec = importlib.util.spec_from_file_location("a223_neu_hif", NEU_PY)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def logistic_theta(h, h_lo, h_hi):
    span = max(h_hi - h_lo, 1e-6)
    h0 = 0.5 * (h_lo + h_hi)
    return 1.0 / (1.0 + np.exp((h - h0) / (span / 6.0)))


def auroc(a, b):
    if min(len(a), len(b)) < 8:
        return np.nan
    y = np.r_[np.ones(len(a)), np.zeros(len(b))]
    return float(roc_auc_score(y, np.r_[a, b]))


def mwu(a, b):
    if min(len(a), len(b)) < 8:
        return np.nan
    return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)


def persist(th):
    return np.mean(th <= T_LOW)


def module_log(S, lib, pos, genes):
    ix = [pos[g] for g in genes if g in pos]
    if not ix:
        return np.zeros(S.shape[0])
    x = np.asarray(S[:, ix].sum(axis=1)).ravel()
    return np.log1p(x * (1e4 / np.maximum(lib, 1.0)))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    neu = load_neu()
    cells = pd.read_csv(CELLS)
    nmask = cells["palak_neutrophil"].eq(True) & cells["spliced_umi"].ge(NEU_UMI)
    tmask = cells["palak_tumor"].eq(True) & cells["spliced_umi"].ge(TUMOR_UMI)
    neu_df = cells[nmask].copy()
    tum_df = cells[tmask].copy()
    dcf_n = neu_df["sample_id"].eq("hypoxia_plus")
    dp_n = neu_df["sample_id"].eq("DP")
    dcf_t = tum_df["sample_id"].eq("hypoxia_plus")
    dp_t = tum_df["sample_id"].eq("DP")
    dn_t = tum_df["sample_id"].eq("DN")

    th_n = neu_df["theta_neu_scale"].to_numpy()
    th_t = tum_df["theta_normoxic"].to_numpy()
    hif_n = neu_df["hif_score_neu_scale"].to_numpy()

    mapping = neu.load_symbols()
    qc = pd.read_csv(GENE_QC)
    use = qc["use_theta"].astype(str).str.lower().isin(["true", "1"])
    genes = qc.loc[use, "gene"].astype(str).tolist()
    w = qc.set_index("gene").loc[genes, "cohens_d_e15_vs_e14"].clip(lower=0).to_numpy(dtype=np.float64)
    w = w / w.sum()
    ref = ad.read_h5ad(E14E15)
    neu.set_symbols(ref, mapping)
    cal = neu.e14_calibration(mapping, genes, w)
    libr = np.asarray(ref.layers["spliced"].sum(axis=1), dtype=np.float64).ravel()
    posr = {str(g): i for i, g in enumerate(ref.var_names.astype(str))}
    log_r = np.log1p(neu.size_normalize(neu.csr_cols(ref.layers["spliced"], [posr[g] for g in genes]), libr))
    h_n_ref = ((log_r - cal["mu_neu"]) / cal["sd_neu"]) @ w
    th_ref = logistic_theta(h_n_ref, cal["h_lo_neu"], cal["h_hi_neu"])
    sample = ref.obs["sample"].astype(str).to_numpy()
    cg = ref.obs["cell_group"].astype(str).to_numpy()
    e14n = (cg == "neutrophil") & (libr >= NEU_UMI) & (sample == "E14S")
    e15n = (cg == "neutrophil") & (libr >= NEU_UMI) & (sample == "E15S")

    a223 = neu.load_merged(mapping)
    scores = neu.annotate(a223)
    Sa = a223.layers["spliced"].tocsr()
    liba = scores["spliced_umi"].to_numpy()
    posa = {str(g): i for i, g in enumerate(a223.var_names.astype(str))}
    neu_a = scores["palak_neutrophil"].to_numpy() & (liba >= NEU_UMI)
    gate = scores["sample_id"].astype(str).to_numpy()
    mpo = module_log(Sa, liba, posa, MPO_EST)
    hifm = module_log(Sa, liba, posa, HIF)
    dcf_a = neu_a & (gate == "hypoxia_plus")
    dp_a = neu_a & (gate == "DP")

    annotated = cells["palak_matched"].eq(True)
    neu_all = cells["palak_neutrophil"].eq(True) & annotated
    dcf_all = cells["sample_id"].eq("hypoxia_plus") & annotated
    dn_all = cells["sample_id"].eq("DN") & annotated
    n_dcf, k_neu_dcf = int(dcf_all.sum()), int((dcf_all & neu_all).sum())
    n_dn, k_neu_dn = int(dn_all.sum()), int((dn_all & neu_all).sum())
    table = np.array([[k_neu_dcf, n_dcf - k_neu_dcf], [k_neu_dn, n_dn - k_neu_dn]])
    or_neu, p_neu = stats.fisher_exact(table)

    tests = []

    def add(name, prediction, observed, support, **extra):
        tests.append({"test": name, "dye_prediction": prediction, "observed": observed, "supports_dye_not_hypoxia": support, **extra})

    p_dcf = persist(th_n[dcf_n.to_numpy()])
    p_dp = persist(th_n[dp_n.to_numpy()])
    p_e14 = persist(th_ref[e14n])
    p_e15 = persist(th_ref[e15n])
    p_tum_dcf = persist(th_t[dcf_t.to_numpy()])
    p_tum_dp = persist(th_t[dp_t.to_numpy()])

    add(
        "DCF+ neutrophils vs never-hypoxic E14 neutrophils (neu-scale θ)",
        "DCF+ ≈ E14 (reverted / HIF-low), not HIF-high",
        f"DCF+ persist {100*p_dcf:.1f}% median θ {np.median(th_n[dcf_n]):.2f}; E14 persist {100*p_e14:.1f}% median θ {np.median(th_ref[e14n]):.2f}; E15 persist {100*p_e15:.1f}%",
        "partial" if p_dcf > p_e14 + 0.1 else "yes",
        auroc_hif_DCF_vs_E14=auroc(hif_n[dcf_n.to_numpy()], h_n_ref[e14n]),
    )
    add(
        "DCF+-only vs DP neutrophils (lactate as orthogonal hypoxia/metabolic mark)",
        "DP more HIF-persistent than DCF+-only if DCF-only is dye",
        f"DCF+ persist {100*p_dcf:.1f}% median θ {np.median(th_n[dcf_n]):.2f}; DP persist {100*p_dp:.1f}% median θ {np.median(th_n[dp_n]):.2f}",
        "yes" if p_dp > p_dcf + 0.1 else "no",
        auroc_hif_DP_vs_DCF=auroc(hif_n[dp_n.to_numpy()], hif_n[dcf_n.to_numpy()]),
        mwu_p=mwu(hif_n[dp_n.to_numpy()], hif_n[dcf_n.to_numpy()]),
    )
    add(
        "Same DCF+ gate: Tumor vs neutrophil HIF-persistent fraction",
        "Tumors HIF-high, neutrophils not, if gate is not a shared hypoxic state",
        f"DCF+ Tumor persist {100*p_tum_dcf:.1f}%; DCF+ neutrophil persist {100*p_dcf:.1f}% (neu-scale). DP Tumor {100*p_tum_dp:.1f}%; DP neutrophil {100*p_dp:.1f}%",
        "yes" if abs(p_tum_dcf - p_dcf) > 0.05 else "no",
    )
    add(
        "Neutrophil overrepresentation in DCF+ vs HIF among DCF+ neutrophils",
        "Lineage enrichment in the FITC gate without a matching HIF-persistent majority",
        f"Neutrophil OR vs DN {or_neu:.2f} (p={p_neu:.1e}); DCF+ neutrophils only {100*p_dcf:.1f}% HIF-persistent",
        "yes" if (or_neu > 2 and p_dcf < 0.5) else "no",
        OR_neu_DCF_vs_DN=or_neu,
        p_fisher=p_neu,
        pct_neu_DCF=100 * k_neu_dcf / n_dcf,
        pct_neu_DN=100 * k_neu_dn / n_dn,
    )
    add(
        "Dye-handling (Mpo/esterase) vs HIF RNA in DCF+ vs DP neutrophils",
        "DCF+-only higher dye/MPO program, not higher HIF",
        f"Mpo/esterase AUROC DCF+>DP {auroc(mpo[dcf_a], mpo[dp_a]):.2f}; HIF module AUROC DCF+>DP {auroc(hifm[dcf_a], hifm[dp_a]):.2f}",
        "yes" if auroc(mpo[dcf_a], mpo[dp_a]) > 0.55 and auroc(hifm[dcf_a], hifm[dp_a]) < 0.5 else "partial",
        auroc_mpo_DCF_vs_DP=auroc(mpo[dcf_a], mpo[dp_a]),
        auroc_hif_DCF_vs_DP=auroc(hifm[dcf_a], hifm[dp_a]),
    )

    verdict = pd.DataFrame(tests)
    yes = (verdict["supports_dye_not_hypoxia"] == "yes").sum()
    partial = (verdict["supports_dye_not_hypoxia"] == "partial").sum()
    verdict.to_csv(OUT / "hypothesis_tests.csv", index=False)
    verdict.to_csv(REPO / "dye_vs_hypoxia_tests.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    data = [
        th_ref[e14n],
        th_ref[e15n],
        th_n[dcf_n.to_numpy()],
        th_n[dp_n.to_numpy()],
    ]
    labels = [
        f"E14 neu\n(n={e14n.sum()})",
        f"E15 neu\n(n={e15n.sum()})",
        f"A223 DCF+\n(n={int(dcf_n.sum())})",
        f"A223 DP\n(n={int(dp_n.sum())})",
    ]
    ax.boxplot(data, labels=labels, showfliers=False)
    ax.axhline(T_LOW, color="#c0392b", ls="--", lw=1)
    ax.axhline(T_HIGH, color="#2980b9", ls="--", lw=1)
    ax.set_ylabel("θ (neutrophil HIF scale; high = reverted)")
    ax.set_title("Dye vs hypoxia: neutrophil θ")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "neu_theta_dcf_vs_hypoxia.png", dpi=180)
    plt.close(fig)

    print(verdict.to_string(index=False))
    print(f"tests supporting dye-not-hypoxia: {yes} yes, {partial} partial / {len(verdict)}")
    del ref, a223


if __name__ == "__main__":
    main()
