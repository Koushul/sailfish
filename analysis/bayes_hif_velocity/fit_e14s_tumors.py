#!/usr/bin/env python3
"""Fit the HIF-α lag model on E14S tumor cells only (single cohort).

No second sample, no neutrophil pool, no forced reverted labels. h and θ are
relative to this cohort's own median. Reads the placed h5ad read-only.
"""
from __future__ import annotations

import argparse
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
from load_h5ad import load_placed
from model import (
    Data,
    fit,
    hypoxia_factor,
    phenotype_calls,
)

DEFAULT_H5AD = "/ix1/ylee/kor11/MC38/E14SE15S/E14SE15S_gex_adt_placed.h5ad"
SAMPLE = "E14S"
LINEAGE = "Tumor"
MIN_SPLICED = 5000.0
ANCHOR_GENES = {"Car9", "Vegfa", "Adm", "Angptl4", "Bnip3", "Ndrg1", "Ddit4"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", default=DEFAULT_H5AD)
    p.add_argument("--out-dir", default=str(HERE / "results"))
    p.add_argument("--n-steps", type=int, default=550)
    p.add_argument("--n-boot", type=int, default=80)
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def gene_mask(U: np.ndarray) -> np.ndarray:
    frac = (U > 0).mean(axis=0)
    return (frac >= 0.05) & (U.mean(axis=0) >= 0.02)


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


def clone_data(data: Data) -> Data:
    return Data(
        U=data.U.copy(),
        S=data.S.copy(),
        L=data.L.copy(),
        theta=data.theta.copy(),
        cycle_s=data.cycle_s.copy(),
        cycle_g2m=data.cycle_g2m.copy(),
        exposed=data.exposed.copy(),
        use_velocity=data.use_velocity.copy(),
        log_kappa0=data.log_kappa0.copy(),
        anchor=data.anchor.copy(),
    )


def build_data(placed) -> tuple[Data, np.ndarray, np.ndarray]:
    keep = (placed.sample == SAMPLE) & (placed.cell_group == LINEAGE) & (placed.L >= MIN_SPLICED)
    n = int(keep.sum())
    if n < 20:
        raise SystemExit(f"too few cells: {n}")
    U = placed.U[keep]
    S = placed.S[keep]
    keep_g = gene_mask(U)
    use_v = velocity_mask(U, S, keep_g)
    cs = placed.cycle_s[keep]
    cg = placed.cycle_g2m[keep]
    cs = cs - float(cs.mean())
    cg = cg - float(cg.mean())
    targets = core_targets()
    anchor = np.array([2.0 if t.mouse in ANCHOR_GENES else 1.0 for t in targets], dtype=np.float64)
    data = Data(
        U=U,
        S=S,
        L=placed.L[keep],
        theta=np.full(n, 0.5),
        cycle_s=cs,
        cycle_g2m=cg,
        exposed=np.zeros(n),
        use_velocity=use_v.astype(np.float64),
        log_kappa0=np.zeros(U.shape[1]),
        anchor=anchor,
    )
    return data, keep, keep_g


def bootstrap_h(data: Data, n_boot: int, seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    n = data.U.shape[0]
    acc = [[] for _ in range(n)]
    frac = np.zeros(n_boot)
    for b in range(n_boot):
        ix = rng.integers(0, n, size=n)
        sub = Data(
            U=data.U[ix],
            S=data.S[ix],
            L=data.L[ix],
            theta=data.theta[ix],
            cycle_s=data.cycle_s[ix],
            cycle_g2m=data.cycle_g2m[ix],
            exposed=data.exposed[ix],
            use_velocity=data.use_velocity,
            log_kappa0=data.log_kappa0,
            anchor=data.anchor,
        )
        fac = hypoxia_factor(sub)
        frac[b] = float(np.mean(fac["h"] >= 1.5))
        for j, i in enumerate(ix):
            acc[i].append(fac["h"][j])
    h_sd = np.array([np.std(a, ddof=1) if len(a) > 1 else np.nan for a in acc])
    p_hi = np.array([np.mean(np.array(a) >= 1.5) if a else np.nan for a in acc])
    return h_sd, p_hi, frac


def entropy(p: np.ndarray) -> np.ndarray:
    q = np.clip(p, 1e-12, 1.0)
    q = q / q.sum(axis=1, keepdims=True)
    return -(q * np.log(q)).sum(axis=1)


def write_md(path: Path, rows: list[str]) -> None:
    path.write_text("\n".join(rows) + "\n")


def main() -> None:
    args = parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    placed = load_placed(args.h5ad)
    data, mask, keep_g = build_data(placed)
    print(
        f"E14S tumors n={data.U.shape[0]} detected={int(keep_g.sum())} "
        f"velocity={int((data.use_velocity > 0.5).sum())} "
        f"Tirosh S={placed.n_s_genes} G2M={placed.n_g2m_genes}",
        flush=True,
    )
    est = fit(clone_data(data), n_steps=args.n_steps, seed=args.seed, use_cycle=True)
    est0 = fit(clone_data(data), n_steps=args.n_steps, seed=args.seed, use_cycle=False)
    pheno = phenotype_calls(est["theta"], force_control_reverted=False)
    pheno0 = phenotype_calls(est0["theta"], force_control_reverted=False)
    idx = np.flatnonzero(mask)
    rho_h_s, _ = spearmanr(est["h"], data.cycle_s)
    rho_h_g, _ = spearmanr(est["h"], data.cycle_g2m)
    rho_xi_s, _ = spearmanr(est["xi_mean"], data.cycle_s)
    rho_xi_g, _ = spearmanr(est["xi_mean"], data.cycle_g2m)
    rho_h0_s, _ = spearmanr(est0["h"], data.cycle_s)
    rho_h_h0, _ = spearmanr(est["h"], est0["h"])
    logL = np.log(np.clip(data.L, 1.0, None))
    rho_h_L, _ = spearmanr(est["h"], logL)

    h_sd, persist_boot, frac_persist_boot = bootstrap_h(data, args.n_boot, args.seed + 7)

    snr = np.abs(est["xi_mean"]) / np.clip(est["xi_sd"], 1e-6, None)
    xi_ci0 = (np.abs(est["xi_mean"]) < 1.96 * est["xi_sd"]).mean()
    H_ph = entropy(est["p_pheno"])
    H_fx = entropy(est["p_flux"])

    qs = pd.Series(pd.qcut(data.cycle_s, 5, duplicates="drop"), dtype="category")
    cycle_rows = []
    for lab in qs.cat.categories:
        m = qs == lab
        cycle_rows.append(
            {
                "s_quintile": str(lab),
                "n": int(m.sum()),
                "mean_h": float(est["h"][m].mean()),
                "frac_persist": float(np.mean(pheno[m] == "persistent")),
                "mean_p_away": float(est["p_away"][m].mean()),
                "mean_cycle_s": float(data.cycle_s[m].mean()),
            }
        )

    cells = pd.DataFrame(
        {
            "cell": placed.cell_id[idx],
            "spliced_umi": data.L,
            "cycle_s": data.cycle_s,
            "cycle_g2m": data.cycle_g2m,
            "h": est["h"],
            "h_sd_boot": h_sd,
            "theta": est["theta"],
            "pheno": pheno,
            "pheno_nocycle": pheno0,
            "h_nocycle": est0["h"],
            "p_persist": est["p_pheno"][:, 0],
            "p_partial": est["p_pheno"][:, 1],
            "p_reverted": est["p_pheno"][:, 2],
            "xi": est["xi_mean"],
            "xi_sd": est["xi_sd"],
            "p_away": est["p_away"],
            "p_toward": est["p_toward"],
            "p_none": est["p_none"],
            "state": est["state"],
            "boot_p_persist": persist_boot,
        }
    )
    genes = pd.DataFrame(
        [
            {
                "mouse": t.mouse,
                "detected": bool(keep_g[i]),
                "use_velocity": bool(data.use_velocity[i] > 0.5),
                "beta": float(est["beta"][i]),
                "lambda": float(est["lambda"][i]),
                "kappa": float(est["kappa"][i]),
                "c_s": float(est["c_s"][i]),
                "c_g2m": float(est["c_g2m"][i]),
            }
            for i, t in enumerate(core_targets())
        ]
    )
    cells.to_csv(out / "e14s_tumor_cells.tsv", sep="\t", index=False)
    genes.to_csv(out / "e14s_tumor_genes.tsv", sep="\t", index=False)
    pd.DataFrame(cycle_rows).to_csv(out / "e14s_tumor_cycle.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(10.2, 8.4))
    ax = axes[0, 0]
    ax.hist(est["h"], bins=25, color="#4d4d4d")
    ax.axvline(1.5, color="#b2182b", ls="--", label="persist h=1.5")
    ax.set_xlabel("h (cohort MAD units)")
    ax.set_ylabel("cells")
    ax.set_title("spliced HIF-target factor")
    ax.legend(frameon=False, fontsize=8)
    ax = axes[0, 1]
    sc = ax.scatter(data.cycle_s, est["h"], c=est["p_away"], s=18, cmap="viridis", vmin=0, vmax=0.4)
    ax.set_xlabel("Tirosh S (cohort-centered)")
    ax.set_ylabel("h")
    ax.set_title(f"h vs cycle S  Spearman={rho_h_s:.2f}")
    fig.colorbar(sc, ax=ax, label="p_away")
    ax = axes[1, 0]
    ax.scatter(est["xi_mean"], est["xi_sd"], s=16, alpha=0.7, c="#2166ac")
    ax.set_xlabel("E[ξ]")
    ax.set_ylabel("Laplace sd(ξ)")
    ax.set_title("lag uncertainty")
    ax.axhline(float(np.median(est["xi_sd"])), color="k", lw=0.6)
    ax = axes[1, 1]
    labs = ["persistent", "partial", "reverted"]
    x = np.arange(2)
    for i, lab in enumerate(labs):
        h = [float(np.mean(pheno == lab)), float(np.mean(pheno0 == lab))]
        ax.bar(x + i * 0.25, h, width=0.24, label=lab)
    ax.set_xticks(x + 0.25)
    ax.set_xticklabels(["with cycle", "no cycle"])
    ax.set_ylabel("fraction")
    ax.set_title("θ-gate vs cycle ablation")
    ax.legend(frameon=False, fontsize=8)
    fig.suptitle("E14S tumors — HIF-α factor (single cohort)", y=0.98)
    fig.tight_layout()
    fig.savefig(out / "e14s_tumor_overview.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    n = len(cells)
    flip = float(np.mean(pheno != pheno0))
    lines = [
        "# E14S tumors (single cohort)",
        "",
        f"Sample `{SAMPLE}`, lineage `{LINEAGE}`, spliced UMI ≥ {MIN_SPLICED:.0f}. "
        "No second library. No forced phenotype. \(h\) is this cohort's own median/MAD. "
        "Persist means ≥1.5 MAD above a typical E14S tumor, not an external never-hypoxic arm.",
        "",
        "## QC",
        "",
        f"- n_cells: {n}",
        f"- median spliced UMI: {float(np.median(data.L)):.0f}",
        f"- detected panel genes: {int(keep_g.sum())} / 32",
        f"- velocity genes: {int((data.use_velocity > 0.5).sum())}",
        "",
        "## Phenotype (θ-gate)",
        "",
        f"- persistent: {float(np.mean(pheno == 'persistent')):.3f} (n={int((pheno == 'persistent').sum())})",
        f"- partial: {float(np.mean(pheno == 'partial')):.3f}",
        f"- reverted: {float(np.mean(pheno == 'reverted')):.3f}",
        f"- GMM persist (argmax): {float(np.mean(est['pheno'] == 'persistent')):.3f}",
        f"- mean p_persist (GMM): {float(est['p_pheno'][:, 0].mean()):.3f}",
        f"- mean phenotype entropy (nats): {float(H_ph.mean()):.3f}",
        "",
        "Relative ranks: about 1.5 MAD above the cohort median is the persist cut. "
        "In a unimodal sample that cut is a tail, not an independent hypoxia class.",
        "",
        "## Bayesian lag uncertainty",
        "",
        f"- σ_v: {float(est['sigma_v']):.3f} (clip is {np.exp(0.5):.3f})",
        f"- median Laplace sd(ξ): {float(np.median(est['xi_sd'])):.3f}",
        f"- mean |ξ|: {float(np.mean(np.abs(est['xi_mean']))):.3f}",
        f"- mean |ξ|/sd: {float(np.mean(snr)):.3f}",
        f"- fraction whose 95% ξ interval includes 0: {xi_ci0:.3f}",
        f"- mean p_none / p_away / p_toward: {float(est['p_none'].mean()):.3f} / {float(est['p_away'].mean()):.3f} / {float(est['p_toward'].mean()):.3f}",
        f"- fraction p_away > 0.5: {float(np.mean(est['p_away'] > 0.5)):.3f}",
        f"- mean flux entropy (nats): {float(H_fx.mean()):.3f}",
        f"- mix π none/away/toward: {est['pi_flux'][0]:.3f} / {est['pi_flux'][1]:.3f} / {est['pi_flux'][2]:.3f}",
        "",
        "Laplace sd(ξ) is the local posterior width of the lag, combined with a spike-slab. "
        "If almost every cell's ξ interval covers 0 and p_none is high, the data do not support directed flux.",
        "",
        "## Bootstrap uncertainty on h",
        "",
        f"- n_boot: {args.n_boot} (cells resampled; ALS factor only)",
        f"- median bootstrap sd(h): {float(np.median(h_sd)):.3f}",
        f"- persist fraction across boots: mean {float(frac_persist_boot.mean()):.3f}, "
        f"2.5–97.5% {float(np.quantile(frac_persist_boot, 0.025)):.3f}–{float(np.quantile(frac_persist_boot, 0.975)):.3f}",
        f"- mean per-cell bootstrap P(h≥1.5): {float(persist_boot.mean()):.3f}",
        "",
        "## Cell cycle",
        "",
        f"- Spearman (h, S): {rho_h_s:.3f}",
        f"- Spearman (h, G2M): {rho_h_g:.3f}",
        f"- Spearman (h, log L): {rho_h_L:.3f}",
        f"- Spearman (ξ, S): {rho_xi_s:.3f}",
        f"- Spearman (ξ, G2M): {rho_xi_g:.3f}",
        f"- Spearman (h_cycle, h_nocycle): {rho_h_h0:.3f}",
        f"- Spearman (h_nocycle, S): {rho_h0_s:.3f}",
        f"- θ-gate label flip with vs without cycle: {flip:.3f}",
        f"- persist with cycle: {float(np.mean(pheno == 'persistent')):.3f}; without: {float(np.mean(pheno0 == 'persistent')):.3f}",
        "",
        "Cycle on spliced z is a gene-specific loading (glycolytic HIF targets track S). "
        "Cycle on unspliced is a gene-specific intercept for lag. ξ is residualized on S/G2M. "
        "Ablation zeros both S and G2M scores.",
        "",
        "### Persist by Tirosh S quintile (with cycle in the model)",
        "",
    ]
    try:
        lines.append(pd.DataFrame(cycle_rows).to_markdown(index=False))
    except Exception:
        lines.append("```\n" + pd.DataFrame(cycle_rows).to_string(index=False) + "\n```")
    top = genes.sort_values("beta", ascending=False).head(8)
    lines += [
        "",
        "## Top factor genes (β)",
        "",
        ", ".join(f"{m}={b:.2f}" for m, b in zip(top.mouse, top.beta)),
        "",
        "## How to read this",
        "",
        "- This is one library of tumor cells. Persist/reverted are tails of that library's HIF-target factor.",
        "- Lag posteriors are Bayesian (Laplace + spike-slab). h bootstrap is frequentist resampling of the spliced factor.",
        "- If cycle ablation barely moves h ranks but persist % moves, the 1.5 MAD cut is sensitive, not the factor.",
        "",
    ]
    write_md(out / "e14s_tumors.md", lines)
    print("wrote", out / "e14s_tumors.md", flush=True)


if __name__ == "__main__":
    main()
