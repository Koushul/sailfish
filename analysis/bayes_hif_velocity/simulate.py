#!/usr/bin/env python3
"""Synthetic spliced/unspliced counts from a HIF-α target panel."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hif_targets import HifTarget, core_targets
from model import Data


def logistic_theta(h: np.ndarray, h_lo: float, h_hi: float) -> np.ndarray:
    span = max(h_hi - h_lo, 1e-6)
    h0 = 0.5 * (h_lo + h_hi)
    tau = span / 6.0
    return 1.0 / (1.0 + np.exp((h - h0) / tau))


@dataclass
class SimConfig:
    n_control: int = 250
    n_exposed: int = 750
    n_core: int | None = None
    n_decoy: int = 0
    dropout: float = 0.12
    cycle_on_hif: float = 0.15
    cycle_on_counts: float = 1.0
    capture_shift: float = 0.5
    lag_scale: float = 1.0
    phi: float = 8.0
    mean_library: float = 8000.0
    scramble_unspliced: bool = False
    hif_spliced_scale: float = 3.0
    cycle_on_spliced: float = 0.18
    phi_spliced: float = 20.0
    seed: int = 0


def _gene_kinetics(genes: list[HifTarget], rng: np.random.Generator) -> dict[str, np.ndarray]:
    g = len(genes)
    w = np.empty(g)
    lam = np.empty(g)
    c_s = np.empty(g)
    expr = np.empty(g)
    drop = np.empty(g)
    kappa = np.empty(g)
    for i, gene in enumerate(genes):
        if gene.role == "glycolysis":
            w[i] = rng.uniform(0.9, 1.4)
            lam[i] = rng.uniform(0.9, 1.5)
            c_s[i] = rng.uniform(0.35, 0.75)
            expr[i] = rng.uniform(0.9, 1.6)
            drop[i] = rng.uniform(0.04, 0.10)
            kappa[i] = np.exp(rng.normal(np.log(0.18), 0.15))
        elif gene.role in ("angiogenesis", "peptide", "pH"):
            w[i] = rng.uniform(1.1, 1.7)
            lam[i] = rng.uniform(0.7, 1.3)
            c_s[i] = rng.uniform(0.05, 0.25)
            expr[i] = rng.uniform(0.25, 0.7)
            drop[i] = rng.uniform(0.15, 0.35)
            kappa[i] = np.exp(rng.normal(np.log(0.10), 0.2))
        elif gene.role == "feedback":
            w[i] = rng.uniform(0.7, 1.1)
            lam[i] = rng.uniform(0.5, 1.0)
            c_s[i] = rng.uniform(0.05, 0.2)
            expr[i] = rng.uniform(0.4, 0.9)
            drop[i] = rng.uniform(0.10, 0.22)
            kappa[i] = np.exp(rng.normal(np.log(0.14), 0.2))
        else:
            w[i] = rng.uniform(0.6, 1.2)
            lam[i] = rng.uniform(0.6, 1.2)
            c_s[i] = rng.uniform(0.10, 0.40)
            expr[i] = rng.uniform(0.45, 1.1)
            drop[i] = rng.uniform(0.08, 0.20)
            kappa[i] = np.exp(rng.normal(np.log(0.14), 0.2))
    w = w / w.sum()
    return {
        "w": w,
        "lambda": lam,
        "c_s": c_s,
        "expr": expr,
        "drop": drop,
        "kappa": kappa,
    }


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
    cycle_s = rng.normal(0.0, 1.0, n)
    cycle_g2m = 0.6 * cycle_s + rng.normal(0.0, 0.8, n)
    cycle_s = cycle_s - cycle_s[exposed < 0.5].mean()
    cycle_g2m = cycle_g2m - cycle_g2m[exposed < 0.5].mean()

    w = np.zeros(n_genes)
    w[:n_core] = kin["w"]
    lam = np.zeros(n_genes)
    lam[:n_core] = kin["lambda"] * cfg.lag_scale
    c_s = np.zeros(n_genes)
    c_s[:n_core] = kin["c_s"] * cfg.cycle_on_counts
    c_g2m = rng.normal(0.0, 0.12, n_genes)
    c_g2m[:n_core] *= cfg.cycle_on_counts
    expr = np.ones(n_genes)
    expr[:n_core] = kin["expr"]
    gene_drop = np.full(n_genes, cfg.dropout)
    gene_drop[:n_core] = np.clip(kin["drop"] * (cfg.dropout / 0.12), 0.0, 0.8)
    kappa = np.exp(rng.normal(np.log(0.08), 0.2, n_genes))
    kappa[:n_core] = kin["kappa"]
    rho0 = np.exp(rng.normal(cfg.capture_shift, 0.15, n_genes))
    phi = np.full(n_genes, cfg.phi)

    if n_decoy:
        c_s[n_core:] = rng.uniform(0.5, 1.0, n_decoy) * cfg.cycle_on_counts
        expr[n_core:] = rng.uniform(0.8, 1.5, n_decoy)
        gene_drop[n_core:] = rng.uniform(0.05, 0.12, n_decoy)
        names = [g.human for g in genes] + [f"DECOY{i+1}" for i in range(n_decoy)]
    else:
        names = [g.human for g in genes]

    h_latent = np.zeros(n)
    h_latent[exposed < 0.5] = rng.normal(-1.55, 0.15, cfg.n_control)
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
        p=[0.22, 0.16, 0.26, 0.12, 0.08, 0.08, 0.08],
    )
    h_exp = np.empty(cfg.n_exposed)
    xi_true = np.zeros(n)
    for i, k in enumerate(kind):
        if k == "persist":
            h_exp[i] = rng.normal(1.7, 0.18)
        elif k == "persistent_exiting":
            h_exp[i] = rng.normal(1.55, 0.18)
            xi_true[cfg.n_control + i] = rng.uniform(0.6, 1.3)
        elif k == "partial":
            h_exp[i] = rng.normal(0.05, 0.12)
        elif k == "transitioning_out":
            h_exp[i] = rng.normal(0.1, 0.12)
            xi_true[cfg.n_control + i] = rng.uniform(0.6, 1.3)
        elif k == "transitioning_in":
            h_exp[i] = rng.normal(0.0, 0.12)
            xi_true[cfg.n_control + i] = -rng.uniform(0.6, 1.3)
        elif k == "reverted":
            h_exp[i] = rng.normal(-1.5, 0.18)
        else:
            h_exp[i] = rng.normal(-1.35, 0.18)
            xi_true[cfg.n_control + i] = -rng.uniform(0.6, 1.3)
    h_latent[exposed > 0.5] = h_exp
    h_latent = h_latent + cfg.cycle_on_hif * cycle_s

    L_bg = rng.lognormal(np.log(cfg.mean_library), 0.35, n)
    L_bg[exposed < 0.5] *= rng.uniform(0.35, 0.55)
    L_bg = np.clip(L_bg, 800, None)

    z = cfg.hif_spliced_scale * h_latent[:, None] * (
        np.clip(w, 1e-8, None) / (w[:n_core].mean() if n_core else 1.0)
    )[None, :]
    z[:, n_core:] = 0.0
    z = z + cfg.cycle_on_spliced * cycle_s[:, None] * (c_s / (np.abs(c_s).mean() + 1e-6))[None, :]
    phi_s = np.full(n_genes, cfg.phi_spliced)
    mean_s = L_bg[:, None] * 0.003 * np.exp(z) * expr[None, :]
    S = rng.negative_binomial(phi_s, np.clip(phi_s / (phi_s + mean_s), 1e-4, 1 - 1e-4))
    L_obs = L_bg + S.sum(axis=1).astype(np.float64)
    L_obs = np.maximum(L_obs, 1.0)

    log1p_s = np.log1p(S * (np.median(L_obs) / L_obs)[:, None])
    mu0 = log1p_s[exposed < 0.5].mean(axis=0)
    sd0 = np.clip(log1p_s[exposed < 0.5].std(axis=0), 1e-6, None)
    zscore = (log1p_s - mu0) / sd0
    h = zscore[:, :n_core] @ w[:n_core] if n_core else zscore.mean(axis=1)
    h_lo = float(np.median(h[exposed < 0.5]))
    q = float(np.quantile(h[exposed > 0.5], 0.75))
    hi = (exposed > 0.5) & (h >= q)
    h_hi = float(np.median(h[hi])) if hi.any() else float(np.median(h[exposed > 0.5]))
    theta = logistic_theta(h, h_lo, h_hi)

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
    U = rng.negative_binomial(phi, np.clip(phi / (phi + mu_u), 1e-4, 1 - 1e-4))
    drop = rng.random(U.shape) < gene_drop[None, :]
    U = np.where(drop, 0, U)
    if cfg.scramble_unspliced:
        for j in range(n_genes):
            U[:, j] = rng.permutation(U[:, j])

    use_velocity = lam > 0.05
    low_th = (exposed > 0.5) & (theta <= np.quantile(theta[exposed > 0.5], 0.2))
    if not np.any(low_th):
        low_th = exposed > 0.5
    kappa_hat = (U[low_th].mean(axis=0) + 1e-3) / (S[low_th].mean(axis=0) + 1e-3)
    log_kappa0 = np.log(np.clip(kappa_hat, 1e-4, 10.0))

    data = Data(
        U=U.astype(np.float64),
        S=S.astype(np.float64),
        L=L_obs,
        theta=theta,
        cycle_s=cycle_s,
        cycle_g2m=cycle_g2m,
        exposed=exposed,
        use_velocity=use_velocity,
        log_kappa0=log_kappa0,
    )
    truth = {
        "xi": xi_true,
        "kind": np.array(["control"] * cfg.n_control + list(kind), dtype=object),
        "kappa": kappa,
        "lambda": lam,
        "rho0": rho0,
        "c_s": c_s,
        "weights": w,
        "theta": theta,
        "h": h,
        "genes": np.array(names, dtype=object),
        "n_core": n_core,
        "n_decoy": n_decoy,
        "cfg": cfg,
    }
    return data, truth
