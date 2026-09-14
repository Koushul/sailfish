#!/usr/bin/env python3
"""Negative-binomial lag model: MAP + Laplace posterior for residual hypoxia lag ξ."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import expit

EPS = 1e-8


def softplus(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    return np.where(x > 20.0, x, np.log1p(np.exp(-np.abs(x))) + np.maximum(x, 0.0))


def d_softplus(x: np.ndarray) -> np.ndarray:
    return expit(x)


@dataclass
class Data:
    U: np.ndarray
    S: np.ndarray
    L: np.ndarray
    theta: np.ndarray
    cycle_s: np.ndarray
    cycle_g2m: np.ndarray
    exposed: np.ndarray
    use_velocity: np.ndarray
    log_kappa0: np.ndarray


def pack_sizes(n: int, g: int) -> list[tuple[str, int]]:
    return [
        ("log_kappa", g),
        ("ell_lambda", g),
        ("c_s", g),
        ("c_g2m", g),
        ("c_r", g),
        ("log_rho0", g),
        ("a_theta", 1),
        ("a_s", 1),
        ("a_g2m", 1),
        ("a_r", 1),
        ("log_sigma_v", 1),
        ("log_phi", g),
        ("xi_hat", n),
    ]


def split(vec: np.ndarray, n: int, g: int) -> dict[str, np.ndarray]:
    out = {}
    i = 0
    for name, size in pack_sizes(n, g):
        out[name] = vec[i : i + size]
        i += size
    return out


def n_params(n: int, g: int) -> int:
    return sum(s for _, s in pack_sizes(n, g))


def init_loc(data: Data) -> np.ndarray:
    n, g = data.U.shape
    p = {
        "log_kappa": data.log_kappa0.copy(),
        "ell_lambda": np.zeros(g),
        "c_s": np.zeros(g),
        "c_g2m": np.zeros(g),
        "c_r": np.zeros(g),
        "log_rho0": np.zeros(g),
        "a_theta": np.zeros(1),
        "a_s": np.zeros(1),
        "a_g2m": np.zeros(1),
        "a_r": np.zeros(1),
        "log_sigma_v": np.array([-2.0]),
        "log_phi": np.full(g, 2.0),
        "xi_hat": np.zeros(n),
    }
    return np.concatenate([p[k] for k, _ in pack_sizes(n, g)])


def xi_from(p: dict, data: Data) -> tuple[np.ndarray, np.ndarray, float]:
    sigma = float(np.exp(np.clip(p["log_sigma_v"][0], -8.0, 2.0)))
    th = data.theta - float(data.theta.mean())
    raw = (
        p["a_theta"][0] * th
        + p["a_s"][0] * data.cycle_s
        + p["a_g2m"][0] * data.cycle_g2m
        + p["a_r"][0] * data.exposed
        + sigma * p["xi_hat"]
    )
    ctrl = data.exposed < 0.5
    cmean = float(raw[ctrl].mean())
    return raw - cmean, raw, sigma


def mu_from(p: dict, data: Data, xi: np.ndarray) -> np.ndarray:
    lam = softplus(p["ell_lambda"]) * data.use_velocity
    kappa = np.exp(np.clip(p["log_kappa"], -8.0, 8.0))
    rho0 = np.exp(np.clip(p["log_rho0"], -4.0, 4.0))
    s_norm = (data.S + 0.5) / data.L[:, None]
    log_rho = np.log(rho0)[None, :] * (1.0 - data.exposed)[:, None]
    log_mu = (
        np.log(data.L)[:, None]
        + log_rho
        + np.log(kappa)[None, :]
        + np.log(np.clip(s_norm, EPS, None))
        + lam[None, :] * xi[:, None]
        + p["c_s"][None, :] * data.cycle_s[:, None]
        + p["c_g2m"][None, :] * data.cycle_g2m[:, None]
        + p["c_r"][None, :] * data.exposed[:, None]
    )
    return np.exp(np.clip(log_mu, -20.0, 20.0)), lam, kappa, rho0


def objective_and_grad(vec: np.ndarray, data: Data) -> tuple[float, np.ndarray]:
    n, g = data.U.shape
    p = split(vec, n, g)
    xi, raw, sigma = xi_from(p, data)
    mu, lam, kappa, rho0 = mu_from(p, data, xi)
    phi = np.exp(np.clip(p["log_phi"], -2.0, 8.0))
    phi_b = phi[None, :]
    resid = (data.U - mu) * phi_b / (phi_b + mu)
    nll = -nb_ll(data.U, mu, phi)
    npri = -log_prior(p, data)
    loss = nll + npri

    dlogmu = resid
    g_xi = (dlogmu * lam[None, :]).sum(axis=1)
    ctrl = data.exposed < 0.5
    n0 = float(ctrl.sum())
    g_xi = g_xi - (g_xi[ctrl].sum() / n0) * ctrl.astype(np.float64)

    th = data.theta - float(data.theta.mean())

    grads = {
        "log_kappa": dlogmu.sum(axis=0) - (p["log_kappa"] - data.log_kappa0) / (0.3**2),
        "ell_lambda": (dlogmu * xi[:, None] * data.use_velocity[None, :]).sum(axis=0)
        * d_softplus(p["ell_lambda"])
        - p["ell_lambda"]
        - 2e6 * p["ell_lambda"] * (~data.use_velocity.astype(bool)),
        "c_s": (dlogmu * data.cycle_s[:, None]).sum(axis=0) - p["c_s"] / (0.3**2),
        "c_g2m": (dlogmu * data.cycle_g2m[:, None]).sum(axis=0) - p["c_g2m"] / (0.3**2),
        "c_r": (dlogmu * data.exposed[:, None]).sum(axis=0) - p["c_r"] / (0.3**2),
        "log_rho0": (dlogmu * (1.0 - data.exposed)[:, None]).sum(axis=0) - p["log_rho0"] / (0.5**2),
        "a_theta": np.array([np.dot(g_xi, th)]) - p["a_theta"] / (0.5**2),
        "a_s": np.array([np.dot(g_xi, data.cycle_s)]) - p["a_s"] / (0.5**2),
        "a_g2m": np.array([np.dot(g_xi, data.cycle_g2m)]) - p["a_g2m"] / (0.5**2),
        "a_r": np.array([np.dot(g_xi, data.exposed)]) - p["a_r"] / (0.5**2),
        "log_sigma_v": np.array([np.dot(g_xi, p["xi_hat"]) * sigma])
        - (p["log_sigma_v"] + 2.0) / (0.5**2),
        "log_phi": d_log_phi(data.U, mu, phi) - (p["log_phi"] - 2.0),
        "xi_hat": g_xi * sigma - p["xi_hat"],
    }

    gvec = np.concatenate([grads[k] for k, _ in pack_sizes(n, g)])
    return float(loss), -gvec


def d_log_phi(k: np.ndarray, mu: np.ndarray, phi: np.ndarray) -> np.ndarray:
    from scipy.special import digamma

    phi_b = phi[None, :]
    term = (
        digamma(k + phi_b)
        - digamma(phi_b)
        + np.log(phi_b)
        - np.log(phi_b + mu)
        + (mu - k) / (phi_b + mu)
    )
    return (term * phi_b).sum(axis=0)


def nb_ll(k: np.ndarray, mu: np.ndarray, phi: np.ndarray) -> float:
    from scipy.special import gammaln

    mu = np.clip(mu, EPS, 1e12)
    phi = np.clip(phi, 1e-3, 1e6)
    phi_b = phi if phi.ndim == 2 else phi[None, :]
    p = phi_b / (phi_b + mu)
    npar = phi_b
    ll = gammaln(k + npar) - gammaln(npar) - gammaln(k + 1.0) + npar * np.log(np.clip(p, EPS, 1.0))
    ll = ll + k * np.log(np.clip(1.0 - p, EPS, 1.0))
    return float(np.sum(ll))


def log_prior(p: dict, data: Data) -> float:
    lp = 0.0
    lp += -0.5 * np.sum(((p["log_kappa"] - data.log_kappa0) / 0.3) ** 2)
    lp += -0.5 * np.sum(p["ell_lambda"] ** 2)
    lp += -0.5 * np.sum((p["c_s"] / 0.3) ** 2 + (p["c_g2m"] / 0.3) ** 2 + (p["c_r"] / 0.3) ** 2)
    lp += -0.5 * np.sum((p["log_rho0"] / 0.5) ** 2)
    for k in ("a_theta", "a_s", "a_g2m", "a_r"):
        lp += -0.5 * (p[k][0] / 0.5) ** 2
    lp += -0.5 * ((p["log_sigma_v"][0] + 2.0) / 0.5) ** 2
    lp += -0.5 * np.sum(((p["log_phi"] - 2.0) / 1.0) ** 2)
    lp += -0.5 * np.sum(p["xi_hat"] ** 2)
    lp += -1e6 * np.sum(p["ell_lambda"][~data.use_velocity.astype(bool)] ** 2)
    return float(lp)


def fit(
    data: Data,
    n_steps: int = 600,
    lr: float = 0.05,
    seed: int = 0,
) -> dict:
    n, g = data.U.shape
    x = init_loc(data)
    m1 = np.zeros_like(x)
    m2 = np.zeros_like(x)
    b1, b2, aeps = 0.9, 0.999, 1e-8
    hist = []
    rng = np.random.default_rng(seed)
    x = x + 0.01 * rng.normal(size=x.shape)
    for t in range(1, n_steps + 1):
        loss, grad = objective_and_grad(x, data)
        m1 = b1 * m1 + (1 - b1) * grad
        m2 = b2 * m2 + (1 - b2) * (grad * grad)
        mh = m1 / (1 - b1**t)
        vh = m2 / (1 - b2**t)
        x = x - lr * mh / (np.sqrt(vh) + aeps)
        hist.append(loss)
    p = split(x, n, g)
    xi, _, sigma = xi_from(p, data)
    mu, lam, _, _ = mu_from(p, data, xi)
    phi = np.exp(np.clip(p["log_phi"], -2.0, 8.0))
    fish = ((phi[None, :] * mu) / (phi[None, :] + mu) * (lam[None, :] ** 2)).sum(axis=1)
    var = 1.0 / np.clip(fish + 1.0, 1e-6, None)
    sd = np.sqrt(var) * sigma
    from scipy.stats import norm

    p_away = 1.0 - norm.cdf(0.0, loc=xi, scale=np.clip(sd, 1e-4, 10.0))
    p_toward = 1.0 - p_away
    return {
        "params": p,
        "xi_mean": xi,
        "p_away": p_away,
        "p_toward": p_toward,
        "xi_sd": sd,
        "loss": hist,
        "lambda": lam,
        "kappa": np.exp(np.clip(p["log_kappa"], -8, 8)),
        "rho_control": np.exp(np.clip(p["log_rho0"], -4, 4)),
        "c_s": p["c_s"],
        "c_g2m": p["c_g2m"],
        "sigma_v": sigma,
        "a_s": float(p["a_s"][0]),
        "a_g2m": float(p["a_g2m"][0]),
    }


def direction_calls(theta: np.ndarray, p_away: np.ndarray, p_toward: np.ndarray, exposed: np.ndarray) -> np.ndarray:
    ctrl = exposed < 0.5
    q_away = float(np.quantile(p_away[ctrl], 0.95))
    q_toward = float(np.quantile(p_toward[ctrl], 0.95))
    out = np.array(["partial"] * theta.size, dtype=object)
    out[theta <= 0.3] = "persistent"
    out[theta >= 0.7] = "reverted"
    away = (p_away >= q_away) & (theta <= 0.7)
    toward = (p_toward >= q_toward) & (theta >= 0.3)
    out[away] = "reverting"
    out[toward] = "inducing"
    return out
