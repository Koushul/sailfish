#!/usr/bin/env python3
"""Realistic synthetic spliced/unspliced counts with known hypoxia lag."""
from __future__ import annotations

import numpy as np

from model import Data


def logistic_theta(h: np.ndarray, h_lo: float, h_hi: float) -> np.ndarray:
    span = max(h_hi - h_lo, 1e-6)
    h0 = 0.5 * (h_lo + h_hi)
    tau = span / 6.0
    return 1.0 / (1.0 + np.exp((h - h0) / tau))


def simulate(
    n_control: int = 250,
    n_exposed: int = 750,
    n_genes: int = 10,
    seed: int = 0,
) -> tuple[Data, dict]:
    rng = np.random.default_rng(seed)
    n = n_control + n_exposed
    exposed = np.r_[np.zeros(n_control), np.ones(n_exposed)]
    cycle_s = rng.normal(0.0, 1.0, n)
    cycle_g2m = 0.6 * cycle_s + rng.normal(0.0, 0.8, n)
    cycle_s = cycle_s - cycle_s[exposed < 0.5].mean()
    cycle_g2m = cycle_g2m - cycle_g2m[exposed < 0.5].mean()

    w = rng.uniform(0.4, 1.5, n_genes)
    w = w / w.sum()
    kappa = np.exp(rng.normal(np.log(0.15), 0.25, n_genes))
    lam = np.r_[rng.uniform(0.8, 1.6, n_genes - 2), 0.0, 0.0]
    use_velocity = lam > 0.05
    c_s = rng.normal(0.0, 0.25, n_genes)
    c_s[:3] = rng.uniform(0.35, 0.7, 3)
    c_g2m = rng.normal(0.0, 0.15, n_genes)
    rho0 = np.exp(rng.normal(0.5, 0.15, n_genes))
    phi = np.full(n_genes, 8.0)

    h_latent = np.zeros(n)
    h_latent[exposed < 0.5] = rng.normal(-1.2, 0.25, n_control)
    kind = rng.choice(
        ["persist", "partial", "reverted", "reverting", "inducing"],
        size=n_exposed,
        p=[0.28, 0.18, 0.34, 0.12, 0.08],
    )
    h_exp = np.empty(n_exposed)
    xi_true = np.zeros(n)
    for i, k in enumerate(kind):
        if k == "persist":
            h_exp[i] = rng.normal(1.4, 0.3)
        elif k == "partial":
            h_exp[i] = rng.normal(0.2, 0.3)
        elif k == "reverted":
            h_exp[i] = rng.normal(-1.1, 0.3)
        elif k == "reverting":
            h_exp[i] = rng.normal(0.4, 0.35)
            xi_true[n_control + i] = rng.uniform(0.6, 1.3)
        else:
            h_exp[i] = rng.normal(-0.2, 0.35)
            xi_true[n_control + i] = -rng.uniform(0.6, 1.3)
    h_latent[exposed > 0.5] = h_exp
    h_latent = h_latent + 0.15 * cycle_s

    L = rng.lognormal(np.log(8000), 0.35, n)
    L[exposed < 0.5] *= rng.uniform(0.35, 0.55)
    L = np.clip(L, 800, None)

    z = h_latent[:, None] * (w / w.mean())[None, :]
    z = z + 0.4 * cycle_s[:, None] * (c_s / (np.abs(c_s).mean() + 1e-6))[None, :]
    mean_s = L[:, None] * np.exp(z - 3.0) * rng.uniform(0.5, 1.5, n_genes)
    S = rng.negative_binomial(phi, np.clip(phi / (phi + mean_s), 1e-4, 1 - 1e-4))
    S = np.maximum(S, 0)
    L_obs = np.maximum(S.sum(axis=1).astype(np.float64), 1.0)

    log1p_s = np.log1p(S * (np.median(L_obs) / L_obs)[:, None])
    mu0 = log1p_s[exposed < 0.5].mean(axis=0)
    sd0 = np.clip(log1p_s[exposed < 0.5].std(axis=0), 1e-6, None)
    h = ((log1p_s - mu0) / sd0) @ w
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
    drop = rng.random(U.shape) < 0.12
    U = np.where(drop, 0, U)

    low_th = (exposed > 0.5) & (theta <= np.quantile(theta[exposed > 0.5], 0.2))
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
        "kind": np.array(["control"] * n_control + list(kind), dtype=object),
        "kappa": kappa,
        "lambda": lam,
        "rho0": rho0,
        "c_s": c_s,
        "weights": w,
        "theta": theta,
        "h": h,
    }
    return data, truth
