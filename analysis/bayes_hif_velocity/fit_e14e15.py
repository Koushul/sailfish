#!/usr/bin/env python3
"""Fit HIF-α lag on E14 (control) vs E15 (exposed), tumors and neutrophils separately.

Reads the placed h5ad read-only. Writes tables and figures under results/.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from hif_targets import core_targets
from load_h5ad import PlacedCounts, load_placed
from model import FLUX_STATES, STATE_NAMES, Data, fit, hypoxia_factor

DEFAULT_H5AD = "/ix1/ylee/kor11/MC38/E14SE15S/E14SE15S_gex_adt_placed.h5ad"
CONTROL_SAMPLE = "E14S"
EXPOSED_SAMPLE = "E15S"
LINEAGES = {
    "Tumor": {"min_spliced": 5000.0},
    "neutrophil": {"min_spliced": 2000.0},
}
ANCHOR_GENES = {"Car9", "Vegfa", "Adm", "Angptl4", "Bnip3", "Ndrg1", "Ddit4"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", default=DEFAULT_H5AD)
    p.add_argument("--out-dir", default=str(Path(__file__).resolve().parent / "results"))
    p.add_argument("--n-steps", type=int, default=550)
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def gene_mask(U: np.ndarray, exposed: np.ndarray) -> np.ndarray:
    frac = (U > 0).mean(axis=0)
    mean_exp = U[exposed > 0.5].mean(axis=0) if np.any(exposed > 0.5) else U.mean(axis=0)
    return (frac >= 0.05) & (mean_exp >= 0.02)


def velocity_mask(U: np.ndarray, S: np.ndarray, keep: np.ndarray) -> np.ndarray:
    use = keep.copy()
    use[S.mean(axis=0) < 0.05] = False
    logk = np.log((U.mean(axis=0) + 1e-3) / (S.mean(axis=0) + 1e-3))
    if keep.sum() >= 5:
        med = np.median(logk[keep])
        mad = np.median(np.abs(logk[keep] - med)) + 1e-6
        use[np.abs(logk - med) > 3.5 * 1.4826 * mad] = False
    use &= keep
    if not np.any(use):
        use = keep.copy()
    return use


def subset_lineage(placed: PlacedCounts, lineage: str, min_spliced: float) -> dict:
    in_lineage = placed.cell_group == lineage
    sample_ok = np.isin(placed.sample, [CONTROL_SAMPLE, EXPOSED_SAMPLE])
    umi_ok = placed.L >= min_spliced
    keep_cells = in_lineage & sample_ok & umi_ok
    exposed = (placed.sample[keep_cells] == EXPOSED_SAMPLE).astype(np.float64)
    U = placed.U[keep_cells]
    S = placed.S[keep_cells]
    keep_g = gene_mask(U, exposed)
    use_v = velocity_mask(U, S, keep_g)
    cs = placed.cycle_s[keep_cells]
    cg = placed.cycle_g2m[keep_cells]
    ctrl = exposed < 0.5
    cs = cs - float(cs[ctrl].mean())
    cg = cg - float(cg[ctrl].mean())
    targets = core_targets()
    anchor = np.array([2.0 if t.mouse in ANCHOR_GENES else 1.0 for t in targets], dtype=np.float64)
    data = Data(
        U=U,
        S=S,
        L=placed.L[keep_cells],
        theta=np.full(int(keep_cells.sum()), 0.5),
        cycle_s=cs,
        cycle_g2m=cg,
        exposed=exposed,
        use_velocity=use_v.astype(np.float64),
        log_kappa0=np.zeros(U.shape[1]),
        anchor=anchor,
    )
    return {
        "data": data,
        "mask": keep_cells,
        "keep_genes": keep_g,
        "n_in_lineage": int(in_lineage.sum()),
        "n_umi_fail": int((in_lineage & sample_ok & ~umi_ok).sum()),
        "n_control": int(ctrl.sum()),
        "n_exposed": int((~ctrl).sum()),
    }


def cell_table(placed: PlacedCounts, mask: np.ndarray, data: Data, est: dict, lineage: str) -> pd.DataFrame:
    idx = np.flatnonzero(mask)
    df = pd.DataFrame(
        {
            "cell": placed.cell_id[idx],
            "lineage": lineage,
            "sample": placed.sample[idx],
            "exposed": (placed.sample[idx] == EXPOSED_SAMPLE).astype(int),
            "spliced_umi": placed.L[idx],
            "h": est["h"],
            "theta": data.theta,
            "xi": est["xi_mean"],
            "xi_sd": est["xi_sd"],
            "p_away": est["p_away"],
            "p_toward": est["p_toward"],
            "p_none": est["p_none"],
            "p_persist": est["p_pheno"][:, 0],
            "p_partial": est["p_pheno"][:, 1],
            "p_reverted": est["p_pheno"][:, 2],
            "pheno": est["pheno"],
            "state": est["state"],
            "historical_state": placed.historical_state[idx],
            "historical_theta": placed.historical_theta[idx],
        }
    )
    for j, name in enumerate(STATE_NAMES):
        df[f"p_{name}"] = est["p_state"][:, j]
    return df


def gene_table(est: dict, keep_g: np.ndarray, use_v: np.ndarray, beta: np.ndarray) -> pd.DataFrame:
    rows = []
    for i, t in enumerate(core_targets()):
        rows.append(
            {
                "mouse": t.mouse,
                "human": t.human,
                "role": t.role,
                "detected": bool(keep_g[i]),
                "use_velocity": bool(use_v[i] > 0.5),
                "lambda": float(est["lambda"][i]),
                "kappa": float(est["kappa"][i]),
                "rho_control": float(est["rho_control"][i]),
                "c_s": float(est["c_s"][i]),
                "c_g2m": float(est["c_g2m"][i]),
                "beta": float(beta[i]),
            }
        )
    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame, lineage: str, n_steps: int, loss0: float, loss1: float, n_genes_v: int) -> pd.DataFrame:
    rows = []
    for sample, sub in df.groupby("sample", sort=False):
        row = {
            "lineage": lineage,
            "sample": sample,
            "n_cells": len(sub),
            "n_velocity_genes": n_genes_v,
            "n_steps": n_steps,
            "loss_start": loss0,
            "loss_final": loss1,
            "mean_h": float(sub["h"].mean()),
            "mean_theta": float(sub["theta"].mean()),
            "mean_xi": float(sub["xi"].mean()),
            "mean_p_away": float(sub["p_away"].mean()),
            "mean_p_toward": float(sub["p_toward"].mean()),
            "frac_hard_flux": float(np.mean(sub["state"].isin(FLUX_STATES))),
        }
        for lab in ("persistent", "partial", "reverted"):
            row[f"frac_pheno_{lab}"] = float(np.mean(sub["pheno"] == lab))
            row[f"mean_p_{lab}"] = float(sub[f"p_{lab}" if lab != "persistent" else "p_persist"].mean())
        for lab in STATE_NAMES:
            row[f"frac_state_{lab}"] = float(np.mean(sub["state"] == lab))
            row[f"mean_p_state_{lab}"] = float(sub[f"p_{lab}"].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def eval_lineage(df: pd.DataFrame, data: Data, est: dict, lineage: str) -> dict:
    ctrl = (df["exposed"] == 0).to_numpy()
    exp = ~ctrl
    cycle_rho, _ = spearmanr(est["xi_mean"], data.cycle_s)
    h_rho_s, _ = spearmanr(est["h"], data.cycle_s)
    out = {
        "lineage": lineage,
        "n_cells": int(len(df)),
        "n_control": int(ctrl.sum()),
        "n_exposed": int(exp.sum()),
        "spearman_xi_cycle_s": float(cycle_rho),
        "spearman_h_cycle_s": float(h_rho_s),
        "mean_abs_xi_control": float(np.mean(np.abs(est["xi_mean"][ctrl]))),
        "mean_p_away_control": float(df.loc[ctrl, "p_away"].mean()),
        "mean_p_toward_control": float(df.loc[ctrl, "p_toward"].mean()),
        "frac_hard_flux_control": float(np.mean(df.loc[ctrl, "state"].isin(FLUX_STATES))),
        "mean_p_away_exposed": float(df.loc[exp, "p_away"].mean()) if exp.any() else float("nan"),
        "mean_p_toward_exposed": float(df.loc[exp, "p_toward"].mean()) if exp.any() else float("nan"),
        "frac_hard_flux_exposed": float(np.mean(df.loc[exp, "state"].isin(FLUX_STATES))) if exp.any() else float("nan"),
        "sigma_v": float(est["sigma_v"]),
        "pi_none": float(est["pi_flux"][0]),
        "pi_away": float(est["pi_flux"][1]),
        "pi_toward": float(est["pi_flux"][2]),
        "mean_lambda_used": float(est["lambda"][data.use_velocity > 0.5].mean())
        if np.any(data.use_velocity > 0.5)
        else float("nan"),
    }
    if exp.any():
        for lab in ("persistent", "partial", "reverted"):
            m = exp & (df["pheno"].to_numpy() == lab)
            out[f"e15_frac_{lab}"] = float(m.mean())
            out[f"e15_n_{lab}"] = int(m.sum())
            if m.any():
                out[f"e15_mean_p_away_{lab}"] = float(df.loc[m, "p_away"].mean())
                out[f"e15_mean_p_toward_{lab}"] = float(df.loc[m, "p_toward"].mean())
    hist = df["historical_theta"].to_numpy()
    tumor_hist = np.isfinite(hist) & exp & (df["lineage"].to_numpy() == "Tumor")
    if tumor_hist.any():
        th = hist[tumor_hist]
        old = np.full(tumor_hist.sum(), "partial", dtype=object)
        old[th <= 0.3] = "persistent"
        old[th >= 0.7] = "reverted"
        new = df.loc[tumor_hist, "pheno"].to_numpy()
        out["hist_n_compared"] = int(tumor_hist.sum())
        out["hist_agree_pheno"] = float(np.mean(old == new))
        for lab in ("persistent", "partial", "reverted"):
            out[f"hist_frac_{lab}"] = float(np.mean(old == lab))
            out[f"new_frac_{lab}_on_hist_cells"] = float(np.mean(new == lab))
    return out


def plot_lineage(df: pd.DataFrame, gene_df: pd.DataFrame, lineage: str, out_dir: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.5))
    samples = [CONTROL_SAMPLE, EXPOSED_SAMPLE]
    pheno_order = ["persistent", "partial", "reverted"]
    colors = {"persistent": "#b2182b", "partial": "#f4a582", "reverted": "#2166ac"}
    ax = axes[0, 0]
    x = np.arange(len(samples))
    bottom = np.zeros(len(samples))
    for lab in pheno_order:
        h = [float(np.mean(df.loc[df["sample"] == s, "pheno"] == lab)) for s in samples]
        ax.bar(x, h, bottom=bottom, color=colors[lab], label=lab)
        bottom += np.array(h)
    ax.set_xticks(x)
    ax.set_xticklabels(["E14 control", "E15 exposed"])
    ax.set_ylabel("fraction")
    ax.set_title(f"{lineage} GMM phenotype")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[0, 1]
    state_show = list(STATE_NAMES)
    cmap = plt.cm.tab20(np.linspace(0, 1, len(state_show)))
    bottom = np.zeros(len(samples))
    for i, lab in enumerate(state_show):
        h = [float(np.mean(df.loc[df["sample"] == s, "state"] == lab)) for s in samples]
        ax.bar(x, h, bottom=bottom, color=cmap[i], label=lab, width=0.7)
        bottom += np.array(h)
    ax.set_xticks(x)
    ax.set_xticklabels(["E14 control", "E15 exposed"])
    ax.set_ylabel("fraction")
    ax.set_title("hard joint state (use posteriors)")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1))

    ax = axes[1, 0]
    for sample, size in [(CONTROL_SAMPLE, 14), (EXPOSED_SAMPLE, 8)]:
        sub = df[df["sample"] == sample]
        ax.scatter(sub["theta"], sub["xi"], s=size, alpha=0.35, label=sample)
    ax.axhline(0.0, color="k", lw=0.6)
    ax.set_xlabel(r"$\theta$ (high = HIF off)")
    ax.set_ylabel(r"$\mathbb{E}[\xi]$ (away $>0$)")
    ax.set_title("spliced phenotype vs lag")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    used = gene_df[gene_df["use_velocity"]]
    order = np.argsort(used["lambda"].to_numpy())
    ax.barh(used["mouse"].to_numpy()[order], used["lambda"].to_numpy()[order], color="#4d4d4d")
    ax.set_xlabel(r"$\lambda_g$")
    ax.set_title("lag genes")
    fig.suptitle(f"E14/E15 HIF-α velocity — {lineage}", y=0.98)
    fig.tight_layout()
    fig.savefig(out_dir / f"e14e15_{lineage.lower()}_overview.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    exp = df[df["exposed"] == 1]
    sc = ax.scatter(exp["p_partial"], exp["p_away"], c=exp["theta"], s=10, cmap="coolwarm_r", vmin=0, vmax=1)
    ax.set_xlabel(r"$p^{\mathrm{partial}}$")
    ax.set_ylabel(r"$p^{\mathrm{away}}$")
    ax.set_title(f"{lineage} E15: partial vs exiting flux")
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(r"$\theta$")
    fig.tight_layout()
    fig.savefig(out_dir / f"e14e15_{lineage.lower()}_partial_flux.png", dpi=140)
    plt.close(fig)


def _fmt_table(df: pd.DataFrame) -> str:
    try:
        return df.to_markdown(index=False)
    except Exception:
        return "```\n" + df.to_string(index=False) + "\n```"


def write_eval_md(path: Path, evals: list[dict], summaries: pd.DataFrame, qc: pd.DataFrame) -> None:
    slim = summaries[
        [
            "lineage",
            "sample",
            "n_cells",
            "frac_pheno_persistent",
            "frac_pheno_partial",
            "frac_pheno_reverted",
            "mean_p_away",
            "mean_p_toward",
            "frac_hard_flux",
        ]
    ]
    lines = [
        "# E14/E15 application",
        "",
        "E14S is the never-hypoxic control; E15S is hypoxia-exposed. Tumors and neutrophils were fit **separately** on the curated HIF-α panel. Spliced library size is cell-wide. Cycle scores are Tirosh S/G2M recomputed on spliced counts and centered on E14. Do not pool lineages. Trust posterior probabilities more than hard `transitioning_out` labels.",
        "",
        "## QC",
        "",
        _fmt_table(qc),
        "",
        "## Sample-level fractions",
        "",
        _fmt_table(slim),
        "",
        "## Diagnostics",
        "",
    ]
    for ev in evals:
        lines.append(f"### {ev['lineage']}")
        lines.append("")
        for k, v in ev.items():
            if k == "lineage":
                continue
            if isinstance(v, float):
                lines.append(f"- `{k}`: {v:.4f}")
            else:
                lines.append(f"- `{k}`: {v}")
        lines.append("")
    lines += [
        "## How to read this",
        "",
        "- Phenotype (`h`, `theta`, GMM persist/partial/reverted) is the spliced HIF-α program. E15 tumors are expected to mix persistent hypoxia with reversion; E14 should sit at the low-HIF end.",
        "- Lag (`xi`, `p_away`, `p_toward`) is residual unspliced after library, capture, and cycle. Synthetic tests recovered rank of `xi` only weakly; report probabilities.",
        "- Neutrophils have lower UMI and a weaker HIF transcriptional program than MC38 tumors. A small persistent fraction there is not evidence they share tumor kinetics.",
        "- Historical `hypoxia_kinetics_state` / `theta_normoxic` on this object were tumor-only gates; they are compared only for tumors and are not used as labels for neutrophils.",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    placed = load_placed(args.h5ad)
    print(
        f"loaded {placed.cell_id.size} cells; Tirosh S genes {placed.n_s_genes}, G2M {placed.n_g2m_genes}",
        flush=True,
    )
    summaries = []
    evals = []
    qc_rows = []
    for lineage, cfg in LINEAGES.items():
        print(f"fitting {lineage} min_spliced={cfg['min_spliced']}", flush=True)
        sub = subset_lineage(placed, lineage, cfg["min_spliced"])
        data = sub["data"]
        qc_rows.append(
            {
                "lineage": lineage,
                "min_spliced_umi": cfg["min_spliced"],
                "n_in_lineage": sub["n_in_lineage"],
                "n_dropped_umi": sub["n_umi_fail"],
                "n_fit": data.U.shape[0],
                "n_control": sub["n_control"],
                "n_exposed": sub["n_exposed"],
                "n_genes": int(data.U.shape[1]),
                "n_detected_genes": int(sub["keep_genes"].sum()),
                "n_velocity_genes": int((data.use_velocity > 0.5).sum()),
                "median_spliced_umi": float(np.median(data.L)),
            }
        )
        est = fit(data, n_steps=args.n_steps, seed=args.seed)
        cells = cell_table(placed, sub["mask"], data, est, lineage)
        genes = gene_table(est, sub["keep_genes"], data.use_velocity, hypoxia_factor(data)["beta"])
        tag = lineage.lower()
        cells.to_csv(out_dir / f"e14e15_{tag}_cells.tsv", sep="\t", index=False)
        genes.to_csv(out_dir / f"e14e15_{tag}_genes.tsv", sep="\t", index=False)
        n_v = int((data.use_velocity > 0.5).sum())
        summ = summarize(cells, lineage, args.n_steps, float(est["loss"][0]), float(est["loss"][-1]), n_v)
        summaries.append(summ)
        ev = eval_lineage(cells, data, est, lineage)
        ev["n_velocity_genes"] = n_v
        ev["loss_final"] = float(est["loss"][-1])
        evals.append(ev)
        plot_lineage(cells, genes, lineage, out_dir)
        e15_pheno = cells.loc[cells["exposed"] == 1, "pheno"].value_counts().to_dict()
        print(
            f"{lineage}: n={len(cells)} E15 pheno {e15_pheno} "
            f"flux_ctrl={ev['frac_hard_flux_control']:.3f} "
            f"xi~cycle={ev['spearman_xi_cycle_s']:.3f}",
            flush=True,
        )
    qc = pd.DataFrame(qc_rows)
    summary = pd.concat(summaries, ignore_index=True)
    ev_df = pd.DataFrame(evals)
    qc.to_csv(out_dir / "e14e15_qc.tsv", sep="\t", index=False)
    summary.to_csv(out_dir / "e14e15_summary.tsv", sep="\t", index=False)
    ev_df.to_csv(out_dir / "e14e15_eval.tsv", sep="\t", index=False)
    write_eval_md(out_dir / "e14e15.md", evals, summary, qc)
    print("wrote", out_dir, flush=True)


if __name__ == "__main__":
    main()
