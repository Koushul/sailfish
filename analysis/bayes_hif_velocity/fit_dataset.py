#!/usr/bin/env python3
"""Fit the HIF-α lag model on a control vs exposed object (one lineage at a time).

Does not write into the h5ad. Skips a lineage if never-hypoxic control n is too small
to identify the control MAD scale.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from hif_targets import core_targets
from load_h5ad import load_counts
from model import (
    FLUX_STATES,
    STATE_NAMES,
    Data,
    fit,
    hypoxia_factor,
    joint_state_probs,
    phenotype_calls,
    state_from_probs,
)

ANCHOR_GENES = {"Car9", "Vegfa", "Adm", "Angptl4", "Bnip3", "Ndrg1", "Ddit4"}
MIN_CONTROL_DEFAULT = 40
RECOMMENDED_CONTROL = 40
PHENO_ORDER = ("persistent", "partial", "reverted")


def _arms(cfg: dict) -> tuple[list[str], list[str]]:
    control = cfg.get("control_samples") or [cfg["control_sample"]]
    exposed = cfg.get("exposed_samples") or [cfg["exposed_sample"]]
    return [str(x) for x in control], [str(x) for x in exposed]


def _barcode16(x) -> str:
    s = str(x)
    if s.startswith("E27_") or s.startswith("E29_"):
        s = s.split("_", 1)[1]
    return s.split("-")[0][:16]


def overlay_qc_tables(placed, tables: list[dict]):
    n = placed.cell_id.size
    kix = {_barcode16(x): i for i, x in enumerate(placed.cell_id)}
    lineage = np.array(["unlabeled"] * n, dtype=object)
    sample = np.array(placed.sample, dtype=object, copy=True)
    hist = np.array(placed.historical_state, dtype=object, copy=True)
    theta = np.array(placed.historical_theta, dtype=np.float64, copy=True)
    n_matched = 0
    for spec in tables:
        df = pd.read_csv(spec["path"])
        lane = spec.get("lane")
        if lane is not None and "lane" in df.columns:
            df = df[df["lane"].astype(str) == str(lane)].copy()
        bc_col = spec.get("barcode_col", "barcode16")
        sample_col = spec.get("sample_col", "sample_id")
        lin_fixed = spec.get("lineage")
        lin_col = spec.get("lineage_col", "cell_group")
        hist_col = spec.get("historical_state_col", "hif_state")
        theta_col = spec.get("historical_theta_col", "theta_normoxic")
        bc = df[bc_col].astype(str).map(_barcode16).to_numpy()
        idx = np.array([kix.get(b, -1) for b in bc], dtype=np.int64)
        ok = idx >= 0
        n_matched += int(ok.sum())
        ii = idx[ok]
        pos = np.flatnonzero(ok)
        if lin_fixed:
            lineage[ii] = lin_fixed
        elif lin_col in df.columns:
            lineage[ii] = df.iloc[pos][lin_col].astype(str).to_numpy()
        if sample_col in df.columns:
            sample[ii] = df.iloc[pos][sample_col].astype(str).to_numpy()
        if hist_col in df.columns:
            hist[ii] = df.iloc[pos][hist_col].astype(str).to_numpy()
        if theta_col in df.columns:
            theta[ii] = pd.to_numeric(df.iloc[pos][theta_col], errors="coerce").to_numpy()
    placed.cell_group = lineage
    placed.sample = sample
    placed.historical_state = hist
    placed.historical_theta = theta
    print(f"qc overlay matched {n_matched} table rows to {int((lineage != 'unlabeled').sum())} cells", flush=True)
    return placed


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True, help="JSON dataset config")
    p.add_argument("--out-dir", default=str(HERE / "results"))
    p.add_argument("--n-steps", type=int, default=550)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--min-control", type=int, default=MIN_CONTROL_DEFAULT)
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


def _onehot_pheno(pheno: np.ndarray) -> np.ndarray:
    return np.column_stack([(pheno == n).astype(np.float64) for n in PHENO_ORDER])


def subset_lineage(
    placed,
    *,
    lineage: str,
    lineage_values: list[str],
    min_spliced: float,
    control_samples: list[str],
    exposed_samples: list[str],
) -> dict:
    in_lineage = np.isin(placed.cell_group, lineage_values)
    sample_ok = np.isin(placed.sample, control_samples + exposed_samples)
    umi_ok = placed.L >= min_spliced
    keep_cells = in_lineage & sample_ok & umi_ok
    exposed = np.isin(placed.sample[keep_cells], exposed_samples).astype(np.float64)
    U = placed.U[keep_cells]
    S = placed.S[keep_cells]
    keep_g = gene_mask(U, exposed)
    use_v = velocity_mask(U, S, keep_g)
    cs = placed.cycle_s[keep_cells]
    cg = placed.cycle_g2m[keep_cells]
    ctrl = exposed < 0.5
    if ctrl.any():
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


def cell_table(placed, mask, data, est, lineage, control_samples, exposed_samples) -> pd.DataFrame:
    idx = np.flatnonzero(mask)
    pheno = phenotype_calls(data.theta, data.exposed, force_control_reverted=True)
    pheno_emp = phenotype_calls(data.theta, data.exposed, force_control_reverted=False)
    p_state = joint_state_probs(_onehot_pheno(pheno), est["p_flux"])
    df = pd.DataFrame(
        {
            "cell": placed.cell_id[idx],
            "lineage": lineage,
            "sample": placed.sample[idx],
            "cell_type": placed.cell_group[idx],
            "exposed": np.isin(placed.sample[idx], exposed_samples).astype(int),
            "spliced_umi": placed.L[idx],
            "h": est["h"],
            "theta": data.theta,
            "xi": est["xi_mean"],
            "xi_sd": est["xi_sd"],
            "p_away": est["p_away"],
            "p_toward": est["p_toward"],
            "p_none": est["p_none"],
            "p_persist_gmm": est["p_pheno"][:, 0],
            "p_partial_gmm": est["p_pheno"][:, 1],
            "p_reverted_gmm": est["p_pheno"][:, 2],
            "p_persist": est["p_pheno"][:, 0],
            "p_partial": est["p_pheno"][:, 1],
            "p_reverted": est["p_pheno"][:, 2],
            "pheno_gmm": est["pheno"],
            "pheno": pheno,
            "pheno_empirical": pheno_emp,
            "state_gmm": est["state"],
            "state": state_from_probs(p_state),
            "historical_state": placed.historical_state[idx],
            "historical_theta": placed.historical_theta[idx],
        }
    )
    for j, name in enumerate(STATE_NAMES):
        df[f"p_{name}"] = p_state[:, j]
    df.attrs["control_samples"] = control_samples
    return df


def gene_table(est, keep_g, use_v, beta) -> pd.DataFrame:
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
        for lab in PHENO_ORDER:
            row[f"frac_pheno_{lab}"] = float(np.mean(sub["pheno"] == lab))
            row[f"frac_pheno_empirical_{lab}"] = float(np.mean(sub["pheno_empirical"] == lab))
            key = "p_persist" if lab == "persistent" else f"p_{lab}"
            row[f"mean_p_{lab}"] = float(sub[key].mean())
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
    logL = np.log(np.clip(data.L, 1.0, None))
    h_rho_L, _ = spearmanr(est["h"], logL)
    out = {
        "lineage": lineage,
        "n_cells": int(len(df)),
        "n_control": int(ctrl.sum()),
        "n_exposed": int(exp.sum()),
        "median_spliced_umi_control": float(np.median(data.L[ctrl])) if ctrl.any() else float("nan"),
        "median_spliced_umi_exposed": float(np.median(data.L[exp])) if exp.any() else float("nan"),
        "umi_ratio_exposed_over_control": (
            float(np.median(data.L[exp]) / max(np.median(data.L[ctrl]), 1.0)) if (exp.any() and ctrl.any()) else float("nan")
        ),
        "spearman_xi_cycle_s": float(cycle_rho),
        "spearman_h_cycle_s": float(h_rho_s),
        "spearman_h_logL": float(h_rho_L),
        "mean_abs_xi_control": float(np.mean(np.abs(est["xi_mean"][ctrl]))) if ctrl.any() else float("nan"),
        "mean_p_away_control": float(df.loc[ctrl, "p_away"].mean()) if ctrl.any() else float("nan"),
        "mean_p_toward_control": float(df.loc[ctrl, "p_toward"].mean()) if ctrl.any() else float("nan"),
        "frac_hard_flux_control": float(np.mean(df.loc[ctrl, "state"].isin(FLUX_STATES))) if ctrl.any() else float("nan"),
        "control_frac_persistent_forced": float(np.mean(df.loc[ctrl, "pheno"] == "persistent")) if ctrl.any() else float("nan"),
        "control_frac_persistent_empirical": float(np.mean(df.loc[ctrl, "pheno_empirical"] == "persistent"))
        if ctrl.any()
        else float("nan"),
        "control_frac_partial_empirical": float(np.mean(df.loc[ctrl, "pheno_empirical"] == "partial")) if ctrl.any() else float("nan"),
        "control_frac_reverted_empirical": float(np.mean(df.loc[ctrl, "pheno_empirical"] == "reverted"))
        if ctrl.any()
        else float("nan"),
        "control_gmm_frac_persistent": float(np.mean(df.loc[ctrl, "pheno_gmm"] == "persistent")) if ctrl.any() else float("nan"),
        "exposed_gmm_frac_persistent": float(np.mean(df.loc[exp, "pheno_gmm"] == "persistent")) if exp.any() else float("nan"),
        "n_control_below_recommended": int(ctrl.sum()) < RECOMMENDED_CONTROL,
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
        for lab in PHENO_ORDER:
            m = exp & (df["pheno"].to_numpy() == lab)
            out[f"exposed_frac_{lab}"] = float(m[exp].mean())
            out[f"exposed_n_{lab}"] = int(m.sum())
            if m.any():
                out[f"exposed_mean_p_away_{lab}"] = float(df.loc[m, "p_away"].mean())
                out[f"exposed_mean_p_toward_{lab}"] = float(df.loc[m, "p_toward"].mean())
    for gate, sub in df.groupby("sample", sort=False):
        out[f"gate_{gate}_n"] = int(len(sub))
        out[f"gate_{gate}_frac_persistent_empirical"] = float(np.mean(sub["pheno_empirical"] == "persistent"))
        out[f"gate_{gate}_mean_h"] = float(sub["h"].mean())
        out[f"gate_{gate}_mean_p_away"] = float(sub["p_away"].mean())
    hist = df["historical_state"].astype(str).to_numpy()
    if exp.any() and np.any(hist != "NA"):
        mapped = np.array(
            [{"control": "reverted", "persistent": "persistent", "partial": "partial", "reverted": "reverted"}.get(x, "other") for x in hist],
            dtype=object,
        )
        use = exp & np.isin(mapped, list(PHENO_ORDER))
        if use.any():
            out["hist_n_compared"] = int(use.sum())
            out["hist_agree_pheno"] = float(np.mean(mapped[use] == df.loc[use, "pheno"].to_numpy()))
            for lab in PHENO_ORDER:
                out[f"hist_frac_{lab}"] = float(np.mean(mapped[use] == lab))
                out[f"new_frac_{lab}_on_hist_cells"] = float(np.mean(df.loc[use, "pheno"].to_numpy() == lab))
    return out


def plot_lineage(df: pd.DataFrame, gene_df: pd.DataFrame, lineage: str, name: str, control_samples: list[str], exposed_samples: list[str], out_dir: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.5))
    samples = [s for s in control_samples + exposed_samples if (df["sample"] == s).any()]
    colors = {"persistent": "#b2182b", "partial": "#f4a582", "reverted": "#2166ac"}
    ax = axes[0, 0]
    x = np.arange(len(samples))
    bottom = np.zeros(len(samples))
    for lab in PHENO_ORDER:
        h = [float(np.mean(df.loc[df["sample"] == s, "pheno"] == lab)) if (df["sample"] == s).any() else 0.0 for s in samples]
        ax.bar(x, h, bottom=bottom, color=colors[lab], label=lab)
        bottom += np.array(h)
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=20, ha="right")
    ax.set_ylabel("fraction")
    ax.set_title(f"{lineage} θ-gate (control forced reverted)")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[0, 1]
    state_show = list(STATE_NAMES)
    cmap = plt.cm.tab20(np.linspace(0, 1, len(state_show)))
    bottom = np.zeros(len(samples))
    for i, lab in enumerate(state_show):
        h = [float(np.mean(df.loc[df["sample"] == s, "state"] == lab)) if (df["sample"] == s).any() else 0.0 for s in samples]
        ax.bar(x, h, bottom=bottom, color=cmap[i], label=lab, width=0.7)
        bottom += np.array(h)
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=20, ha="right")
    ax.set_ylabel("fraction")
    ax.set_title("hard joint state (use posteriors)")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1))

    ax = axes[1, 0]
    for sample in samples:
        sub = df[df["sample"] == sample]
        size = 14 if sample in control_samples else 8
        if len(sub):
            ax.scatter(sub["theta"], sub["xi"], s=size, alpha=0.35, label=sample)
    ax.axhline(0.0, color="k", lw=0.6)
    ax.set_xlabel(r"$\theta$ (high = HIF off)")
    ax.set_ylabel(r"$\mathbb{E}[\xi]$ (away $>0$)")
    ax.set_title("spliced phenotype vs lag")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    used = gene_df[gene_df["use_velocity"]]
    if len(used):
        order = np.argsort(used["lambda"].to_numpy())
        ax.barh(used["mouse"].to_numpy()[order], used["lambda"].to_numpy()[order], color="#4d4d4d")
    ax.set_xlabel(r"$\lambda_g$")
    ax.set_title("lag genes")
    fig.suptitle(f"{name} HIF-α velocity — {lineage}", y=0.98)
    fig.tight_layout()
    fig.savefig(out_dir / f"{name}_{lineage.lower()}_overview.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def _fmt_table(df: pd.DataFrame) -> str:
    try:
        return df.to_markdown(index=False)
    except Exception:
        return "```\n" + df.to_string(index=False) + "\n```"


def write_report(path: Path, name: str, cfg: dict, evals: list[dict], summaries: pd.DataFrame, qc: pd.DataFrame, skips: list[dict]) -> None:
    cols = [
        "lineage",
        "sample",
        "n_cells",
        "frac_pheno_persistent",
        "frac_pheno_partial",
        "frac_pheno_reverted",
        "frac_pheno_empirical_persistent",
        "mean_p_away",
        "mean_p_toward",
        "frac_hard_flux",
    ]
    slim = summaries[cols] if len(summaries) and all(c in summaries.columns for c in cols) else summaries
    lines = [
        f"# {name} application",
        "",
        cfg.get("description", ""),
        "",
        f"Control `{cfg.get('control_samples') or cfg.get('control_sample')}`; exposed `{cfg.get('exposed_samples') or cfg.get('exposed_sample')}`. Forced control labels are design constraints. Empirical θ-gates on control are the diagnostic for circular persist calls. DN n below {RECOMMENDED_CONTROL} means the control MAD for \(h\) is poorly identified — report empirical persist and do not over-read exposed persist.",
        "",
        "## QC",
        "",
        _fmt_table(qc) if len(qc) else "_no lineages fitted_",
        "",
        "## Sample-level fractions",
        "",
        _fmt_table(slim) if len(slim) else "_none_",
        "",
        "## Skipped",
        "",
    ]
    if skips:
        lines.append(_fmt_table(pd.DataFrame(skips)))
    else:
        lines.append("None.")
    lines += ["", "## Diagnostics", ""]
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
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    cfg = json.loads(Path(args.config).read_text())
    name = cfg["name"]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    placed = load_counts(
        cfg["h5ad"],
        sample_col=cfg.get("sample_col", "sample"),
        lineage_col=cfg.get("lineage_col"),
        gene_symbol_col=cfg.get("gene_symbol_col"),
        historical_state_col=cfg.get("historical_state_col"),
        historical_theta_col=cfg.get("historical_theta_col"),
        require_all_panel=bool(cfg.get("require_all_panel", False)),
    )
    if cfg.get("qc_tables"):
        overlay_qc_tables(placed, cfg["qc_tables"])
    print(
        f"{name}: {placed.cell_id.size} cells; Tirosh S {placed.n_s_genes} G2M {placed.n_g2m_genes}; "
        f"missing panel {placed.missing_panel}",
        flush=True,
    )
    control_samples, exposed_samples = _arms(cfg)
    min_control = int(cfg.get("min_control", args.min_control))
    summaries = []
    evals = []
    qc_rows = []
    skips = []
    for lineage, lcfg in cfg["lineages"].items():
        values = lcfg.get("values", [lineage])
        min_spliced = float(lcfg["min_spliced"])
        print(f"subset {lineage} values={values} min_spliced={min_spliced}", flush=True)
        sub = subset_lineage(
            placed,
            lineage=lineage,
            lineage_values=values,
            min_spliced=min_spliced,
            control_samples=control_samples,
            exposed_samples=exposed_samples,
        )
        n_ctrl = sub["n_control"]
        n_exp = sub["n_exposed"]
        row = {
            "lineage": lineage,
            "min_spliced_umi": min_spliced,
            "n_in_lineage": sub["n_in_lineage"],
            "n_dropped_umi": sub["n_umi_fail"],
            "n_fit": sub["data"].U.shape[0],
            "n_control": n_ctrl,
            "n_exposed": n_exp,
            "n_genes": int(sub["data"].U.shape[1]),
            "n_detected_genes": int(sub["keep_genes"].sum()),
            "n_velocity_genes": int((sub["data"].use_velocity > 0.5).sum()),
            "median_spliced_umi": float(np.median(sub["data"].L)) if sub["data"].U.shape[0] else float("nan"),
            "status": "fit",
            "n_control_below_recommended": n_ctrl < RECOMMENDED_CONTROL,
        }
        if n_ctrl < min_control:
            row["status"] = "skip_n_control"
            qc_rows.append(row)
            skips.append(
                {
                    "lineage": lineage,
                    "reason": "n_control below minimum for control MAD / xi pin",
                    "n_control": n_ctrl,
                    "n_exposed": n_exp,
                    "min_control": min_control,
                }
            )
            print(f"SKIP {lineage}: n_control={n_ctrl} < {min_control}", flush=True)
            continue
        if n_exp < 20:
            row["status"] = "skip_n_exposed"
            qc_rows.append(row)
            skips.append(
                {
                    "lineage": lineage,
                    "reason": "n_exposed too small",
                    "n_control": n_ctrl,
                    "n_exposed": n_exp,
                    "min_control": min_control,
                }
            )
            print(f"SKIP {lineage}: n_exposed={n_exp}", flush=True)
            continue
        if n_ctrl < RECOMMENDED_CONTROL:
            print(
                f"WARN {lineage}: n_control={n_ctrl} < {RECOMMENDED_CONTROL}; "
                "control MAD / xi pin is poorly identified",
                flush=True,
            )
        qc_rows.append(row)
        data = sub["data"]
        est = fit(data, n_steps=args.n_steps, seed=args.seed)
        cells = cell_table(placed, sub["mask"], data, est, lineage, control_samples, exposed_samples)
        genes = gene_table(est, sub["keep_genes"], data.use_velocity, hypoxia_factor(data)["beta"])
        tag = lineage.lower().replace(" ", "_")
        cells.to_csv(out_dir / f"{name}_{tag}_cells.tsv", sep="\t", index=False)
        genes.to_csv(out_dir / f"{name}_{tag}_genes.tsv", sep="\t", index=False)
        n_v = int((data.use_velocity > 0.5).sum())
        summ = summarize(cells, lineage, args.n_steps, float(est["loss"][0]), float(est["loss"][-1]), n_v)
        summaries.append(summ)
        ev = eval_lineage(cells, data, est, lineage)
        ev["n_velocity_genes"] = n_v
        ev["loss_final"] = float(est["loss"][-1])
        evals.append(ev)
        plot_lineage(cells, genes, lineage, name, control_samples, exposed_samples, out_dir)
        print(
            f"{lineage}: n={len(cells)} emp_ctrl_persist={ev['control_frac_persistent_empirical']:.3f} "
            f"exp_persist={ev.get('exposed_frac_persistent', float('nan')):.3f} "
            f"p_away ctrl/exp={ev['mean_p_away_control']:.3f}/{ev['mean_p_away_exposed']:.3f}",
            flush=True,
        )
    qc = pd.DataFrame(qc_rows)
    summary = pd.concat(summaries, ignore_index=True) if summaries else pd.DataFrame()
    ev_df = pd.DataFrame(evals) if evals else pd.DataFrame()
    qc.to_csv(out_dir / f"{name}_qc.tsv", sep="\t", index=False)
    if len(summary):
        summary.to_csv(out_dir / f"{name}_summary.tsv", sep="\t", index=False)
    if len(ev_df):
        ev_df.to_csv(out_dir / f"{name}_eval.tsv", sep="\t", index=False)
    write_report(out_dir / f"{name}.md", name, cfg, evals, summary, qc, skips)
    print("wrote", out_dir / f"{name}.md", flush=True)


if __name__ == "__main__":
    main()
