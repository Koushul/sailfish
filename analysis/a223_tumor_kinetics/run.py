#!/usr/bin/env python3
"""Tumor 1-D HIF-down axis on A223 E27/E29 (θ-only persist vs revert)."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

THETA_LO = 0.3
THETA_HI = 0.7
TUMOR_UMI_MIN = 5000
NEU_DIR = Path("/ix1/ylee/shared/sailfish/analysis/a223_neutrophil_hypoxia")

GATE_ORDER = ["DN", "hypoxia_plus", "DP", "lactate_plus"]
GATE_LABEL = {
    "DN": "DN",
    "hypoxia_plus": "DCF+",
    "DP": "DP",
    "lactate_plus": "lactate+",
}


def theta_only_state(theta: np.ndarray) -> np.ndarray:
    out = np.full(theta.shape, "partial", dtype=object)
    out[theta <= THETA_LO] = "persistent"
    out[theta >= THETA_HI] = "reverted"
    return out


def pct_table(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    n = df.groupby(group_cols, observed=True).size().rename("n_qc")
    counts = (
        df.groupby(group_cols + ["theta_state"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    for col in ("persistent", "partial", "reverted"):
        if col not in counts.columns:
            counts[col] = 0
    counts = counts[["persistent", "partial", "reverted"]]
    out = counts.join(n)
    med = df.groupby(group_cols, observed=True)["theta_normoxic"].median().rename("median_theta")
    out = out.join(med)
    for col in ("persistent", "partial", "reverted"):
        out[f"{col}_pct"] = 100.0 * out[col] / out["n_qc"]
    return out.reset_index()


def plot_stacked(df: pd.DataFrame, title: str, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    x = np.arange(len(df))
    persist = df["persistent_pct"].to_numpy()
    partial = df["partial_pct"].to_numpy()
    reverted = df["reverted_pct"].to_numpy()
    ax.bar(x, persist, color="#c0392b", label="persistent (θ≤0.3)")
    ax.bar(x, partial, bottom=persist, color="#f4d03f", label="partial")
    ax.bar(x, reverted, bottom=persist + partial, color="#2980b9", label="reverted (θ≥0.7)")
    ax.set_xticks(x)
    ax.set_xticklabels(df["label"].tolist(), rotation=25, ha="right")
    ax.set_ylabel("% of QC Tumor (spliced UMI ≥ 5000)")
    ax.set_ylim(0, 100)
    ax.set_title(title)
    ax.legend(frameon=False, loc="upper right", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_theta_hist(df: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), sharey=True)
    for ax, lane in zip(axes, ["E27", "E29"]):
        sub = df[df["sample"] == lane]
        ax.hist(sub["theta_normoxic"], bins=40, color="#34495e", alpha=0.85)
        ax.axvline(THETA_LO, color="#c0392b", ls="--", lw=1)
        ax.axvline(THETA_HI, color="#2980b9", ls="--", lw=1)
        ax.set_title(f"{lane} Tumor  (n={len(sub)})")
        ax.set_xlabel("θ_normoxic (tumor HIF-down scale)")
    axes[0].set_ylabel("cells")
    fig.suptitle("A223 Tumor 1-D axis (high θ = reverted / low HIF)", y=1.02)
    fig.tight_layout()
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--cells",
        default="/ix1/ylee/kor11/A223/neutrophil_hypoxia/cells.csv",
    )
    p.add_argument(
        "--e14e15-placed",
        default="/ix1/ylee/kor11/MC38/E14SE15S/E14SE15S_gex_adt_placed.h5ad",
    )
    p.add_argument(
        "--gene-qc",
        default="/ix1/ylee/shared/sailfish/analysis/a223_neutrophil_hypoxia/tumor_theta_gene_qc.csv",
    )
    p.add_argument("--out-dir", default="/ix1/ylee/kor11/A223/tumor_kinetics")
    args = p.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    cells = pd.read_csv(args.cells)
    cells["sample"] = cells["lane"]
    cells["ocm_gate"] = cells["sample_id"]
    cells["umi_spliced"] = cells["spliced_umi"]
    tumor = cells[
        (cells["palak_broad"] == "Tumor")
        & (cells["umi_spliced"] >= TUMOR_UMI_MIN)
        & cells["theta_normoxic"].notna()
    ].copy()
    tumor["theta_state"] = theta_only_state(tumor["theta_normoxic"].to_numpy())

    by_lane = pct_table(tumor, ["sample"])
    by_gate = pct_table(tumor, ["ocm_gate"])
    by_lane_gate = pct_table(tumor, ["sample", "ocm_gate"])
    by_subtype = pct_table(tumor, ["palak_cell_type"])
    pooled = pct_table(tumor.assign(all="all"), ["all"])

    by_lane.to_csv(out / "by_lane.csv", index=False)
    by_gate.to_csv(out / "by_gate.csv", index=False)
    by_lane_gate.to_csv(out / "by_lane_gate.csv", index=False)
    by_subtype.to_csv(out / "by_palak_subtype.csv", index=False)
    pooled.to_csv(out / "pooled.csv", index=False)
    tumor.to_csv(out / "tumor_qc.csv", index=False)

    gene_qc = pd.read_csv(args.gene_qc)
    gene_qc.to_csv(out / "tumor_theta_gene_qc.csv", index=False)

    import importlib.util
    import anndata as ad

    spec = importlib.util.spec_from_file_location("a223_neu_hif", NEU_DIR / "run.py")
    neu = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(neu)

    mapping = neu.load_symbols()
    use = gene_qc["use_theta"].astype(str).str.lower().isin(["true", "1"])
    hif_genes = gene_qc.loc[use, "gene"].astype(str).tolist()
    w = (
        gene_qc.set_index("gene")
        .loc[hif_genes, "cohens_d_e15_vs_e14"]
        .clip(lower=0)
        .to_numpy(dtype=np.float64)
    )
    w = w / w.sum()
    cal = neu.e14_calibration(mapping, hif_genes, w)

    ref = ad.read_h5ad(args.e14e15_placed)
    neu.set_symbols(ref, mapping)
    lib = np.asarray(ref.layers["spliced"].sum(axis=1), dtype=np.float64).ravel()
    var_pos = {str(g): i for i, g in enumerate(ref.var_names.astype(str))}
    cols = [var_pos[g] for g in hif_genes]
    log_s = np.log1p(neu.size_normalize(neu.csr_cols(ref.layers["spliced"], cols), lib))
    h = ((log_s - cal["mu"]) / cal["sd"]) @ w
    theta_e = neu.logistic_theta(h, cal["h_lo"], cal["h_hi"])
    e15_mask = (
        (ref.obs["sample"].astype(str) == "E15S")
        & (ref.obs["cell_group"].astype(str) == "Tumor")
        & (lib >= TUMOR_UMI_MIN)
    )
    e15t = pd.DataFrame({"theta_normoxic": theta_e[np.asarray(e15_mask)]})
    e15t["theta_state"] = theta_only_state(e15t["theta_normoxic"].to_numpy())
    del ref
    e15_row = {
        "cohort": "E15S_MC38",
        "n_qc": int(len(e15t)),
        "persistent": int((e15t["theta_state"] == "persistent").sum()),
        "partial": int((e15t["theta_state"] == "partial").sum()),
        "reverted": int((e15t["theta_state"] == "reverted").sum()),
        "median_theta": float(e15t["theta_normoxic"].median()),
    }
    for k in ("persistent", "partial", "reverted"):
        e15_row[f"{k}_pct"] = 100.0 * e15_row[k] / e15_row["n_qc"]

    a223_row = {
        "cohort": "A223_E27E29",
        "n_qc": int(len(tumor)),
        "persistent": int((tumor["theta_state"] == "persistent").sum()),
        "partial": int((tumor["theta_state"] == "partial").sum()),
        "reverted": int((tumor["theta_state"] == "reverted").sum()),
        "median_theta": float(tumor["theta_normoxic"].median()),
    }
    for k in ("persistent", "partial", "reverted"):
        a223_row[f"{k}_pct"] = 100.0 * a223_row[k] / a223_row["n_qc"]

    compare = pd.DataFrame([e15_row, a223_row])
    compare.to_csv(out / "compare_e15_theta_only.csv", index=False)

    lane_plot = by_lane.copy()
    lane_plot["label"] = lane_plot["sample"]
    plot_stacked(
        lane_plot,
        "A223 Tumor: persistent vs reverted (1-D HIF-down axis)",
        out / "stacked_by_lane.png",
    )
    gate_plot = by_gate.copy()
    gate_plot["ocm_gate"] = pd.Categorical(gate_plot["ocm_gate"], GATE_ORDER)
    gate_plot = gate_plot.sort_values("ocm_gate")
    gate_plot["label"] = gate_plot["ocm_gate"].map(GATE_LABEL)
    plot_stacked(
        gate_plot,
        "A223 Tumor 1-D states by OCM gate (pooled E27+E29)",
        out / "stacked_by_gate.png",
    )
    plot_theta_hist(tumor, out / "theta_hist_by_lane.png")

    comp_plot = pd.DataFrame(
        [
            {
                "label": "E15S MC38\n(θ-only)",
                "persistent_pct": e15_row["persistent_pct"],
                "partial_pct": e15_row["partial_pct"],
                "reverted_pct": e15_row["reverted_pct"],
            },
            {
                "label": "A223 E27+E29\n(θ-only)",
                "persistent_pct": a223_row["persistent_pct"],
                "partial_pct": a223_row["partial_pct"],
                "reverted_pct": a223_row["reverted_pct"],
            },
        ]
    )
    plot_stacked(
        comp_plot,
        "Tumor 1-D persist vs revert: E15S vs A223",
        out / "stacked_vs_e15.png",
    )
    sub_plot = by_subtype.sort_values("n_qc", ascending=False).copy()
    sub_plot["label"] = sub_plot["palak_cell_type"].str.replace("Tumor ", "", regex=False)
    plot_stacked(
        sub_plot,
        "A223 Tumor 1-D states by Tumor subtype",
        out / "stacked_by_subtype.png",
    )

    repo_results = Path("/ix1/ylee/shared/sailfish/analysis/a223_tumor_kinetics/results")
    repo_results.mkdir(parents=True, exist_ok=True)
    for name, df in (
        ("by_lane.csv", by_lane),
        ("by_gate.csv", by_gate),
        ("by_lane_gate.csv", by_lane_gate),
        ("by_palak_subtype.csv", by_subtype),
        ("pooled.csv", pooled),
        ("compare_e15_theta_only.csv", compare),
    ):
        df.to_csv(repo_results / name, index=False)

    lines = [
        "# A223 tumor 1-D hypoxia persistence",
        "",
        "Same Tumor HIF-down axis as E14/E15 (`θ_normoxic`), **θ-only** gates",
        f"(persistent ≤ {THETA_LO}, reverted ≥ {THETA_HI}). No velocity.",
        "QC: annotated `Tumor`, spliced UMI ≥ 5000. High θ = low HIF / reverted.",
        "",
        f"- A223 pooled QC Tumor **n={a223_row['n_qc']}**: "
        f"**{a223_row['persistent_pct']:.1f}% persistent**, "
        f"{a223_row['partial_pct']:.1f}% partial, "
        f"**{a223_row['reverted_pct']:.1f}% reverted**.",
        f"- E15S MC38 Tumor (same θ-only gates, n={e15_row['n_qc']}): "
        f"{e15_row['persistent_pct']:.1f}% / {e15_row['partial_pct']:.1f}% / "
        f"{e15_row['reverted_pct']:.1f}%.",
        "",
        "E15 kinetic labels (velocity) put some θ-reverted cells in inducing/",
        "memory; θ-only is the fair comparison here.",
        "",
        "Do not use `hif_state` from the neutrophil run (neu-scale).",
        "Use `theta_normoxic` + these gates.",
    ]
    (out / "README.md").write_text("\n".join(lines) + "\n")
    print(compare.to_string(index=False))
    print(by_lane.to_string(index=False))
    print(by_gate.to_string(index=False))
    print(by_subtype.to_string(index=False))


if __name__ == "__main__":
    main()
