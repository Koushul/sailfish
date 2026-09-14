#!/usr/bin/env python3
"""Synthetic spliced/unspliced counts from a HIF-α target panel."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hif_targets import HifTarget, core_targets
from model import Data


@dataclass
class SimConfig:
    n_control: int = 280
    n_exposed: int = 920
    n_core: int | None = None
    n_decoy: int = 4
    n_silent: int = 1
    dropout: float = 0.14
    cycle_on_hif: float = 0.2
    cycle_on_counts: float = 1.0
    capture_shift: float = 0.4
    lag_scale: float = 1.0
    phi: float = 7.5
    mean_library: float = 7000.0
    scramble_unspliced: bool = False
    hif_spliced_scale: float = 2.4
    cycle_on_spliced: float = 0.24
    phi_spliced: float = 14.0
    ambient: float = 0.015
    overlap: float = 0.32
    seed: int = 0


def _gene_kinetics(genes: list[HifTarget], rng: np.random.Generator) -> dict[str, np.ndarray]:
    g = len(genes)
    w = np.empty(g)
    lam = np.empty(g)
    c_s = np.empty(g)
    expr = np.empty(g)
    drop = np.empty(g)
    kappa = np.empty(g)
    anchor = np.empty(g)
    for i, gene in enumerate(genes):
        if gene.role == "glycolysis":
            w[i] = rng.uniform(0.9, 1.4)
            lam[i] = rng.uniform(0.9, 1.55)
            c_s[i] = rng.uniform(0.35, 0.85)
            expr[i] = rng.uniform(0.9, 1.7)
            drop[i] = rng.uniform(0.05, 0.12)
            kappa[i] = np.exp(rng.normal(np.log(0.16), 0.2))
            anchor[i] = 0.85
        elif gene.role in ("angiogenesis", "peptide", "pH"):
            w[i] = rng.uniform(1.1, 1.8)
            lam[i] = rng.uniform(0.8, 1.5)
            c_s[i] = rng.uniform(0.04, 0.22)
            expr[i] = rng.uniform(0.2, 0.65)
            drop[i] = rng.uniform(0.22, 0.45)
            kappa[i] = np.exp(rng.normal(np.log(0.09), 0.25))
            anchor[i] = 1.9
        elif gene.role == "feedback":
            w[i] = rng.uniform(0.7, 1.2)
            lam[i] = rng.uniform(0.5, 1.1)
            c_s[i] = rng.uniform(0.05, 0.22)
            expr[i] = rng.uniform(0.35, 0.9)
            drop[i] = rng.uniform(0.14, 0.28)
            kappa[i] = np.exp(rng.normal(np.log(0.13), 0.2))
            anchor[i] = 1.3
        else:
            w[i] = rng.uniform(0.6, 1.25)
            lam[i] = rng.uniform(0.55, 1.25)
            c_s[i] = rng.uniform(0.08, 0.4)
            expr[i] = rng.uniform(0.4, 1.1)
            drop[i] = rng.uniform(0.12, 0.28)
            kappa[i] = np.exp(rng.normal(np.log(0.13), 0.22))
            anchor[i] = 1.35
    w = w / w.sum()
    return {
        "w": w,
        "lambda": lam,
        "c_s": c_s,
        "expr": expr,
        "drop": drop,
        "kappa": kappa,
        "anchor": anchor,
    }


def _lag(rng: np.random.Generator, scale: float) -> float:
    mag = rng.uniform(0.28, 0.5) if rng.random() < 0.3 else rng.uniform(0.7, 1.25)
    return mag * scale


def simulate(cfg: SimConfig | None = None, **kwargs) -> tuple[Data, dict]:
    if cfg is None:
        cfg = SimConfig(**kwargs)
    elif kwargs:
        raise TypeError("pass either SimConfig or keyword args, not both")
    rng = np.random.default_rng(cfg.seed)
    genes = core_targets()
    if cfg.n_core is not None:
        genes = genes[: max(3, min(cfg.n_core, len(genes)))]
    kin = _gene_kinetics(genes, rng)
    n_core = len(genes)
    n_decoy = int(cfg.n_decoy)
    n_genes = n_core + n_decoy
    n = cfg.n_control + cfg.n_exposed
    exposed = np.r_[np.zeros(cfg.n_control), np.ones(cfg.n_exposed)]

    phase = rng.choice([0, 1, 2], size=n, p=[0.55, 0.25, 0.20])
    cycle_s = (phase == 1).astype(np.float64) + 0.22 * (phase == 2) + rng.normal(0.0, 0.32, n)
    cycle_g2m = (phase == 2).astype(np.float64) + 0.18 * (phase == 1) + rng.normal(0.0, 0.32, n)
    cycle_s = cycle_s - cycle_s[exposed < 0.5].mean()
    cycle_g2m = cycle_g2m - cycle_g2m[exposed < 0.5].mean()

    w = np.zeros(n_genes)
    w[:n_core] = kin["w"]
    lam = np.zeros(n_genes)
    lam[:n_core] = kin["lambda"] * cfg.lag_scale
    c_s = np.zeros(n_genes)
    c_s[:n_core] = kin["c_s"] * cfg.cycle_on_counts
    c_g2m = rng.normal(0.0, 0.14, n_genes)
    c_g2m[:n_core] *= cfg.cycle_on_counts
    expr = np.ones(n_genes)
    expr[:n_core] = kin["expr"]
    gene_drop = np.full(n_genes, cfg.dropout)
    gene_drop[:n_core] = np.clip(kin["drop"] * (cfg.dropout / 0.14), 0.0, 0.55)
    kappa = np.exp(rng.normal(np.log(0.08), 0.25, n_genes))
    kappa[:n_core] = kin["kappa"]
    rho0 = np.exp(rng.normal(cfg.capture_shift, 0.2, n_genes))
    phi_u = np.clip(rng.lognormal(np.log(cfg.phi), 0.25, n_genes), 2.0, 20.0)
    anchor = np.full(n_genes, 0.2)
    anchor[:n_core] = kin["anchor"]

    if cfg.n_silent and n_core > 6:
        silent = rng.choice(n_core, size=min(cfg.n_silent, n_core // 3), replace=False)
        w[silent] = 0.0
        lam[silent] = 0.0
        if w[:n_core].sum() > 0:
            w[:n_core] = w[:n_core] / w[:n_core].sum()

    if n_decoy:
        c_s[n_core:] = rng.uniform(0.55, 1.15, n_decoy) * cfg.cycle_on_counts
        expr[n_core:] = rng.uniform(0.7, 1.6, n_decoy)
        gene_drop[n_core:] = rng.uniform(0.06, 0.16, n_decoy)
        names = [g.human for g in genes] + [f"DECOY{i+1}" for i in range(n_decoy)]
    else:
        names = [g.human for g in genes]

    sd = cfg.overlap
    kind = rng.choice(
        [
            "persist",
            "partial",
            "reverted",
            "transitioning_out",
            "transitioning_in",
            "persistent_exiting",
            "reverted_entering",
        ],
        size=cfg.n_exposed,
        p=[0.24, 0.14, 0.36, 0.09, 0.07, 0.05, 0.05],
    )
    h_latent = np.zeros(n)
    h_latent[exposed < 0.5] = rng.normal(-1.25, sd * 0.7, cfg.n_control)
    h_exp = np.empty(cfg.n_exposed)
    xi_true = np.zeros(n)
    for i, k in enumerate(kind):
        if k == "persist":
            h_exp[i] = rng.normal(1.2, sd)
        elif k == "persistent_exiting":
            h_exp[i] = rng.normal(1.05, sd)
            xi_true[cfg.n_control + i] = _lag(rng, cfg.lag_scale)
        elif k == "partial":
            h_exp[i] = rng.normal(0.08, sd * 0.9)
        elif k == "transitioning_out":
            h_exp[i] = rng.normal(0.12, sd * 0.9)
            xi_true[cfg.n_control + i] = _lag(rng, cfg.lag_scale)
        elif k == "transitioning_in":
            h_exp[i] = rng.normal(0.0, sd * 0.9)
            xi_true[cfg.n_control + i] = -_lag(rng, cfg.lag_scale)
        elif k == "reverted":
            h_exp[i] = rng.normal(-1.2, sd)
        else:
            h_exp[i] = rng.normal(-1.05, sd)
            xi_true[cfg.n_control + i] = -_lag(rng, cfg.lag_scale)
    h_latent[exposed > 0.5] = h_exp
    h_latent = h_latent + cfg.cycle_on_hif * cycle_s

    L_bg = rng.lognormal(np.log(cfg.mean_library), 0.5, n)
    L_bg[exposed < 0.5] *= rng.uniform(0.4, 0.7)
    L_bg = np.clip(L_bg, 400, None)

    z = cfg.hif_spliced_scale * h_latent[:, None] * (
        np.clip(w, 1e-8, None) / (w[:n_core].mean() if n_core else 1.0)
    )[None, :]
    z[:, n_core:] = 0.0
    z = z + cfg.cycle_on_spliced * cycle_s[:, None] * (c_s / (np.abs(c_s).mean() + 1e-6))[None, :]
    phi_s = np.full(n_genes, cfg.phi_spliced)
    mean_s = L_bg[:, None] * 0.0024 * np.exp(z) * expr[None, :]
    S = rng.negative_binomial(phi_s, np.clip(phi_s / (phi_s + mean_s), 1e-4, 1 - 1e-4))
    L_obs = np.maximum(L_bg + S.sum(axis=1).astype(np.float64), 1.0)

    s_norm = (S + 0.5) / L_obs[:, None]
    log_mu_u = (
        np.log(L_obs)[:, None]
        + np.log(rho0)[None, :] * (1.0 - exposed)[:, None]
        + np.log(kappa)[None, :]
        + np.log(s_norm)
        + lam[None, :] * xi_true[:, None]
        + c_s[None, :] * cycle_s[:, None]
        + c_g2m[None, :] * cycle_g2m[:, None]
    )
    mu_u = np.exp(np.clip(log_mu_u, -20, 12))
    U = rng.negative_binomial(phi_u, np.clip(phi_u / (phi_u + mu_u), 1e-4, 1 - 1e-4))
    amb = rng.poisson(np.clip(cfg.ambient * mu_u.mean(axis=0), 0.0, 8.0), size=U.shape)
    U = U + amb
    drop = rng.random(U.shape) < gene_drop[None, :]
    U = np.where(drop, 0, U)
    if cfg.scramble_unspliced:
        for j in range(n_genes):
            U[:, j] = rng.permutation(U[:, j])

    use_velocity = lam > 0.05
    data = Data(
        U=U.astype(np.float64),
        S=S.astype(np.float64),
        L=L_obs,
        theta=np.full(n, 0.5),
        cycle_s=cycle_s,
        cycle_g2m=cycle_g2m,
        exposed=exposed,
        use_velocity=use_velocity,
        log_kappa0=np.log(np.clip(kappa, 1e-4, 10.0)),
        anchor=anchor,
    )
    truth = {
        "xi": xi_true,
        "kind": np.array(["control"] * cfg.n_control + list(kind), dtype=object),
        "kappa": kappa,
        "lambda": lam,
        "rho0": rho0,
        "c_s": c_s,
        "weights": w,
        "genes": np.array(names, dtype=object),
        "n_core": n_core,
        "n_decoy": n_decoy,
        "cfg": cfg,
        "phase": phase,
    }
    return data, truth
