#!/usr/bin/env python3
"""HIF-α program factor + spike-slab unspliced lag (MAP, Laplace, mixture posterior)."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.special import expit, logsumexp
from scipy.stats import norm

EPS = 1e-8
STATE_NAMES = (
    "persistent",
    "persistent_exiting",
    "persistent_deepening",
    "partial",
    "transitioning_out",
    "transitioning_in",
    "reverted",
    "reverted_entering",
)
FLUX_STATES = (
    "persistent_exiting",
    "persistent_deepening",
    "transitioning_out",
    "transitioning_in",
    "reverted_entering",
)


def softplus(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    return np.where(x > 20.0, x, np.log1p(np.exp(-np.abs(x))) + np.maximum(x, 0.0))


def d_softplus(x: np.ndarray) -> np.ndarray:
    return expit(x)


def softmax(logp: np.ndarray, axis: int = -1) -> np.ndarray:
    a = logp - np.max(logp, axis=axis, keepdims=True)
    e = np.exp(a)
    return e / np.clip(e.sum(axis=axis, keepdims=True), EPS, None)


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
    anchor: np.ndarray = field(default_factory=lambda: np.array([]))


def pack_sizes(n: int, g: int) -> list[tuple[str, int]]:
    return [
        ("log_kappa", g),
        ("ell_lambda", g),
        ("c_s", g),
        ("c_g2m", g),
        ("c_r", g),
        ("log_rho0", g),
        ("log_sigma_v", 1),
        ("log_phi", g),
        ("ell_omega", g),
        ("flux_logits", 3),
        ("xi_hat", n),
    ]


def split(vec: np.ndarray, n: int, g: int) -> dict[str, np.ndarray]:
    out = {}
    i = 0
    for name, size in pack_sizes(n, g):
        out[name] = vec[i : i + size]
        i += size
    return out


def init_loc(data: Data) -> np.ndarray:
    n, g = data.U.shape
    p = {
        "log_kappa": data.log_kappa0.copy(),
        "ell_lambda": np.full(g, 0.4),
        "c_s": np.zeros(g),
        "c_g2m": np.zeros(g),
        "c_r": np.zeros(g),
        "log_rho0": np.zeros(g),
        "log_sigma_v": np.array([-0.2]),
        "log_phi": np.full(g, 2.0),
        "ell_omega": np.full(g, -1.2),
        "flux_logits": np.array([0.8, 0.0, 0.0]),
        "xi_hat": np.zeros(n),
    }
    return np.concatenate([p[k] for k, _ in pack_sizes(n, g)])


def cycle_design(data: Data) -> np.ndarray:
    x = np.column_stack([data.cycle_s, data.cycle_g2m])
    return x - x.mean(axis=0, keepdims=True)


def has_contrast(exposed: np.ndarray) -> bool:
    e = np.asarray(exposed)
    return bool(np.any(e < 0.5) and np.any(e > 0.5))


def residualize_cycle(v: np.ndarray, data: Data) -> np.ndarray:
    x = cycle_design(data)
    if float(np.max(np.abs(x))) < 1e-12:
        return v
    coef, _, _, _ = np.linalg.lstsq(x, v, rcond=None)
    return v - x @ coef


def mix_params(p: dict) -> tuple[np.ndarray, float, float, float]:
    pi = softmax(p["flux_logits"][None, :])[0]
    return pi, 0.7, 0.28, 0.55


def xi_from(p: dict, data: Data) -> tuple[np.ndarray, float]:
    sigma = float(np.exp(np.clip(p["log_sigma_v"][0], -3.0, 0.5)))
    raw = residualize_cycle(sigma * p["xi_hat"], data)
    if has_contrast(data.exposed):
        ctrl = data.exposed < 0.5
        return raw - float(raw[ctrl].mean()), sigma
    return raw - float(raw.mean()), sigma


def mixture_logp(xi: np.ndarray, pi: np.ndarray, mu: float, t0: float, t1: float):
    means = np.array([0.0, mu, -mu])
    var = np.array([t0**2, t1**2, t1**2])
    x = xi[:, None]
    log_comp = np.log(pi + EPS) - 0.5 * np.log(2.0 * np.pi * var) - 0.5 * (x - means) ** 2 / var
    log_mix = logsumexp(log_comp, axis=1)
    resp = softmax(log_comp, axis=1)
    grad = (resp * (means - x) / var).sum(axis=1)
    return log_mix, grad, resp


def mu_from(p: dict, data: Data, xi: np.ndarray):
    lam = softplus(p["ell_lambda"]) * data.use_velocity
    kappa = np.exp(np.clip(p["log_kappa"], -8.0, 8.0))
    rho0 = np.exp(np.clip(p["log_rho0"], -4.0, 4.0))
    s_norm = (data.S + 0.5) / data.L[:, None]
    if has_contrast(data.exposed):
        log_rho = np.log(rho0)[None, :] * (1.0 - data.exposed)[:, None]
    else:
        log_rho = 0.0
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


def spliced_z(data: Data) -> np.ndarray:
    cpm = 1e4 * data.S / np.clip(data.L, 1.0, None)[:, None]
    log1p_s = np.log1p(cpm)
    ref = data.exposed < 0.5 if has_contrast(data.exposed) else np.ones(data.exposed.size, dtype=bool)
    mu0 = log1p_s[ref].mean(axis=0)
    sd_std = log1p_s[ref].std(axis=0)
    mad = 1.4826 * np.median(np.abs(log1p_s[ref] - mu0), axis=0)
    sd0 = np.clip(np.maximum(sd_std, mad), 0.2, None)
    return np.clip((log1p_s - mu0) / sd0, -6.0, 6.0)


def library_covariate(data: Data) -> np.ndarray:
    ll = np.log(np.clip(data.L, 1.0, None))
    if not has_contrast(data.exposed):
        return (ll - ll.mean()) / (float(ll.std()) + 1e-6)
    ctrl = data.exposed < 0.5
    out = np.empty(ll.shape[0], dtype=np.float64)
    out[ctrl] = ll[ctrl] - ll[ctrl].mean()
    expm = ~ctrl
    if np.any(expm):
        out[expm] = ll[expm] - ll[expm].mean()
    den = float(np.std(ll[ctrl])) + 1e-6
    return out / den


def _balanced_cell_weights(exposed: np.ndarray) -> np.ndarray:
    if not has_contrast(exposed):
        return np.ones(exposed.size, dtype=np.float64)
    ctrl = exposed < 0.5
    n = exposed.size
    n0 = max(int(ctrl.sum()), 1)
    n1 = max(int((~ctrl).sum()), 1)
    w = np.empty(n, dtype=np.float64)
    w[ctrl] = 0.5 * n / n0
    w[~ctrl] = 0.5 * n / n1
    return w


def _logit(p: float) -> float:
    p = float(np.clip(p, 1e-6, 1.0 - 1e-6))
    return float(np.log(p / (1.0 - p)))


def theta_from_control_h(h: np.ndarray, exposed: np.ndarray) -> dict:
    """Map MAD-scaled h to θ. h=0 is the reference-cell median.

    θ=0.88 at h=0; persist (θ=0.3) is 1.5 MADs above that typical cell.
    In a two-sample fit the reference is the r=0 arm; in a single cohort
    it is this cohort's own median (relative ranks, not an external control).
    """
    del exposed
    h_on = 1.5
    t_off, t_on = 0.88, 0.3
    tau = max(h_on / (_logit(t_off) - _logit(t_on)), 0.15)
    h0 = tau * _logit(t_off)
    theta = expit(-(h - h0) / tau)
    return {"theta": theta, "h_lo": 0.0, "h_hi": h_on, "h0": h0, "tau": tau}


def hypoxia_factor(data: Data, n_iter: int = 50, use_cycle: bool = True) -> dict:
    z = spliced_z(data)
    n, g = z.shape
    cs = data.cycle_s if use_cycle else np.zeros(n)
    cg = data.cycle_g2m if use_cycle else np.zeros(n)
    corr = np.array(
        [abs(np.corrcoef(z[:, j], cs)[0, 1]) if z[:, j].std() > 1e-8 else 1.0 for j in range(g)]
    )
    anchor = data.anchor if data.anchor.size == g else np.ones(g)
    w = anchor / (0.3 + corr)
    w = w / np.clip(w.sum(), EPS, None)
    h = z @ w
    h = (h - h.mean()) / (h.std() + EPS)
    ones = np.ones(n)
    lib = library_covariate(data)
    sw = np.sqrt(_balanced_cell_weights(data.exposed))[:, None]
    for _ in range(n_iter):
        H = np.column_stack([h, cs, cg, lib, ones])
        coef = np.linalg.lstsq(H * sw, z * sw, rcond=None)[0]
        beta = np.maximum(coef[0], 0.0)
        resid = (
            z
            - cs[:, None] * coef[1]
            - cg[:, None] * coef[2]
            - lib[:, None] * coef[3]
            - coef[4]
        )
        denom = float(np.dot(beta, beta) + 1e-6)
        h = resid @ beta / denom
        h = h - h.mean()
        scale = h.std()
        if scale > 1e-6:
            h = h / scale
            beta = beta * scale
    beta = np.clip(beta, 0.0, 8.0)
    ctrl = data.exposed < 0.5
    if has_contrast(data.exposed) and h[ctrl].mean() > h[~ctrl].mean():
        h = -h
    h_als = h.copy()
    ref = ctrl if has_contrast(data.exposed) else np.ones(n, dtype=bool)
    h = h - float(np.median(h[ref]))
    mad0 = float(np.median(np.abs(h[ref]))) + 1e-6
    h = h / (1.4826 * mad0)
    h = np.clip(h, -8.0, 8.0)
    mapped = theta_from_control_h(h, data.exposed)
    ll = np.log(np.clip(data.L, 1.0, None))
    mapped["corr_h_logL"] = float(np.corrcoef(h, ll)[0, 1]) if h.size > 2 else 0.0
    mapped["corr_h_cycle_s"] = float(np.corrcoef(h, data.cycle_s)[0, 1]) if h.size > 2 else 0.0
    mapped["lib"] = lib
    mapped["mad0"] = mad0
    return {"h": h, "h_als": h_als, "beta": beta, **mapped}


def phenotype_probs(h: np.ndarray, exposed: np.ndarray, n_iter: int = 40) -> np.ndarray:
    contrast = has_contrast(exposed)
    ctrl = exposed < 0.5
    if contrast:
        mu = np.array(
            [float(np.median(h[ctrl])), float(np.median(h)), float(np.quantile(h[~ctrl], 0.85))]
        )
    else:
        mu = np.quantile(h, [0.2, 0.5, 0.8]).astype(np.float64)
    mu.sort()
    var = np.full(3, max(float(h.var()), 0.05))
    pi = np.array([0.4, 0.3, 0.3])
    boost = np.array([8.0, 1.0, 0.4])
    for _ in range(n_iter):
        log_comp = np.log(pi + EPS) - 0.5 * np.log(2 * np.pi * var) - 0.5 * (h[:, None] - mu) ** 2 / var
        if contrast:
            log_comp[ctrl] += np.log(boost)
        r = softmax(log_comp, axis=1)
        nk = np.clip(r.sum(axis=0), 1e-3, None)
        pi = nk / nk.sum()
        mu = (r * h[:, None]).sum(axis=0) / nk
        var = np.clip((r * (h[:, None] - mu) ** 2).sum(axis=0) / nk, 0.02, 4.0)
        order = np.argsort(mu)
        mu, var, pi = mu[order], var[order], pi[order]
    log_comp = np.log(pi + EPS) - 0.5 * np.log(2 * np.pi * var) - 0.5 * (h[:, None] - mu) ** 2 / var
    resp = softmax(log_comp, axis=1)
    return np.column_stack([resp[:, 2], resp[:, 1], resp[:, 0]])


def annotate_phenotype(data: Data, use_cycle: bool = True) -> dict:
    fac = hypoxia_factor(data, use_cycle=use_cycle)
    data.theta = fac["theta"]
    p_pheno = phenotype_probs(fac["h_als"], data.exposed)
    return {**fac, "p_pheno": p_pheno}


def objective_and_grad(vec: np.ndarray, data: Data) -> tuple[float, np.ndarray]:
    n, g = data.U.shape
    p = split(vec, n, g)
    xi, sigma = xi_from(p, data)
    mu, lam, _, _ = mu_from(p, data, xi)
    phi = np.exp(np.clip(p["log_phi"], -2.0, 8.0))
    omega = expit(p["ell_omega"])
    nll, dlogmu, d_om = zinb_nll_and_dlogmu(data.U, mu, phi, omega)
    npri = -log_prior(p, data)
    loss = nll + npri

    g_xi = (dlogmu * lam[None, :]).sum(axis=1)
    if has_contrast(data.exposed):
        ctrl = data.exposed < 0.5
        n0 = float(ctrl.sum())
        g_xi = g_xi - (g_xi[ctrl].sum() / n0) * ctrl.astype(np.float64)
    else:
        g_xi = g_xi - g_xi.mean()
    g_xi = residualize_cycle(g_xi, data)

    if has_contrast(data.exposed):
        g_c_r = (dlogmu * data.exposed[:, None]).sum(axis=0) - p["c_r"] / (0.3**2)
        g_rho = (dlogmu * (1.0 - data.exposed)[:, None]).sum(axis=0) - p["log_rho0"] / (0.5**2)
    else:
        g_c_r = -p["c_r"] / (0.3**2)
        g_rho = -p["log_rho0"] / (0.5**2)

    grads = {
        "log_kappa": dlogmu.sum(axis=0) - (p["log_kappa"] - data.log_kappa0) / (0.3**2),
        "ell_lambda": (dlogmu * xi[:, None] * data.use_velocity[None, :]).sum(axis=0)
        * d_softplus(p["ell_lambda"])
        - p["ell_lambda"]
        - 2e6 * p["ell_lambda"] * (~data.use_velocity.astype(bool)),
        "c_s": (dlogmu * data.cycle_s[:, None]).sum(axis=0) - p["c_s"] / (0.3**2),
        "c_g2m": (dlogmu * data.cycle_g2m[:, None]).sum(axis=0) - p["c_g2m"] / (0.3**2),
        "c_r": g_c_r,
        "log_rho0": g_rho,
        "log_sigma_v": np.array([np.dot(g_xi, p["xi_hat"]) * sigma]) - (p["log_sigma_v"] + 0.2) / (0.7**2),
        "log_phi": d_log_phi(data.U, mu, phi) - (p["log_phi"] - 2.0),
        "ell_omega": d_om * omega * (1.0 - omega) - (p["ell_omega"] + 1.2),
        "flux_logits": -(p["flux_logits"] - np.array([0.8, 0.0, 0.0])) / 1.5,
        "xi_hat": g_xi * sigma - p["xi_hat"],
    }
    gvec = np.concatenate([grads[k] for k, _ in pack_sizes(n, g)])
    return float(loss), -gvec


def zinb_nll_and_dlogmu(U: np.ndarray, mu: np.ndarray, phi: np.ndarray, omega: np.ndarray):
    phi_b = phi[None, :]
    om = np.clip(omega, 1e-5, 1 - 1e-5)
    om_b = om[None, :]
    p0 = np.exp(phi_b * (np.log(phi_b) - np.log(phi_b + mu)))
    p0 = np.clip(p0, EPS, 1.0)
    zero = U <= 0.0
    pos = ~zero
    mix0 = om_b + (1.0 - om_b) * p0
    ll = np.where(pos, np.log(1.0 - om_b), np.log(np.clip(mix0, EPS, None)))
    nll = -float(np.sum(ll) + nb_ll_pos(U, mu, phi, pos))
    dlog = np.zeros_like(mu)
    resid = (U - mu) * phi_b / (phi_b + mu)
    dlog = np.where(pos, resid, 0.0)
    w0 = ((1.0 - om_b) * p0) / np.clip(mix0, EPS, None)
    dlogp0 = -phi_b * mu / (phi_b + mu)
    dlog = np.where(zero, w0 * dlogp0, dlog)
    d_om = ((1.0 - p0) / np.clip(mix0, EPS, None) * zero).sum(axis=0)
    d_om += -(pos.sum(axis=0) / np.clip(1.0 - om, EPS, None))
    return nll, dlog, d_om


def nb_ll_pos(k: np.ndarray, mu: np.ndarray, phi: np.ndarray, pos: np.ndarray) -> float:
    from scipy.special import gammaln

    mu = np.clip(mu, EPS, 1e12)
    phi = np.clip(phi, 1e-3, 1e6)
    phi_b = phi[None, :]
    p = phi_b / (phi_b + mu)
    npar = phi_b
    ll = gammaln(k + npar) - gammaln(npar) - gammaln(k + 1.0) + npar * np.log(np.clip(p, EPS, 1.0))
    ll = ll + k * np.log(np.clip(1.0 - p, EPS, 1.0))
    return float(np.sum(ll[pos]))


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
    lp += -0.5 * ((p["log_sigma_v"][0] + 0.2) / 0.7) ** 2
    lp += -0.5 * np.sum(((p["log_phi"] - 2.0) / 1.0) ** 2)
    lp += -0.5 * np.sum(p["xi_hat"] ** 2)
    lp += -0.5 * np.sum(((p["ell_omega"] + 1.2) / 1.0) ** 2)
    lp += -1e6 * np.sum(p["ell_lambda"][~data.use_velocity.astype(bool)] ** 2)
    return float(lp)


def flux_posterior(xi: np.ndarray, sd: np.ndarray, pi: np.ndarray, mu: float, t0: float, t1: float):
    means = np.array([0.0, mu, -mu])
    var = np.array([t0**2, t1**2, t1**2]) + sd[:, None] ** 2
    logp = np.log(pi + EPS) - 0.5 * np.log(2 * np.pi * var) - 0.5 * (xi[:, None] - means) ** 2 / var
    return softmax(logp, axis=1)


def joint_state_probs(p_pheno: np.ndarray, p_flux: np.ndarray) -> np.ndarray:
    p_per, p_par, p_rev = p_pheno[:, 0], p_pheno[:, 1], p_pheno[:, 2]
    p_none, p_away, p_tow = p_flux[:, 0], p_flux[:, 1], p_flux[:, 2]
    cols = [
        p_per * p_none,
        p_per * p_away,
        p_per * p_tow,
        p_par * p_none,
        p_par * p_away,
        p_par * p_tow,
        p_rev * (p_none + p_away),
        p_rev * p_tow,
    ]
    p = np.column_stack(cols)
    p = p / np.clip(p.sum(axis=1, keepdims=True), EPS, None)
    return p


def state_from_probs(p_state: np.ndarray) -> np.ndarray:
    idx = np.argmax(p_state, axis=1)
    return np.array([STATE_NAMES[i] for i in idx], dtype=object)


def phenotype_from_probs(p_pheno: np.ndarray) -> np.ndarray:
    names = np.array(["persistent", "partial", "reverted"], dtype=object)
    return names[np.argmax(p_pheno, axis=1)]


def fit(data: Data, n_steps: int = 550, lr: float = 0.04, seed: int = 0, use_cycle: bool = True) -> dict:
    if not use_cycle:
        data.cycle_s = np.zeros_like(data.cycle_s)
        data.cycle_g2m = np.zeros_like(data.cycle_g2m)
    ph = annotate_phenotype(data, use_cycle=True)
    if has_contrast(data.exposed):
        ctrl = data.exposed < 0.5
        low_th = (~ctrl) & (data.theta <= np.quantile(data.theta[~ctrl], 0.2))
        if not np.any(low_th):
            low_th = ~ctrl
    else:
        low_th = data.theta <= np.quantile(data.theta, 0.2)
    kappa_hat = (data.U[low_th].mean(axis=0) + 1e-3) / (data.S[low_th].mean(axis=0) + 1e-3)
    data.log_kappa0 = np.log(np.clip(kappa_hat, 1e-4, 10.0))

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
    xi, sigma = xi_from(p, data)
    mu, lam, _, _ = mu_from(p, data, xi)
    phi = np.exp(np.clip(p["log_phi"], -2.0, 8.0))
    fish = ((phi[None, :] * mu) / (phi[None, :] + mu) * (lam[None, :] ** 2)).sum(axis=1)
    var = 1.0 / np.clip(fish + 1.0 / max(sigma**2, 1e-4), 1e-6, None)
    sd = np.sqrt(var)
    pi = np.array([0.7, 0.15, 0.15])
    mu_f, t0, t1 = 0.7, 0.28, 0.55
    for _ in range(25):
        p_flux = flux_posterior(xi, sd, pi, mu_f, t0, t1)
        pi = 0.85 * p_flux.mean(axis=0) + 0.15 * np.array([0.7, 0.15, 0.15])
        pi = pi / pi.sum()
    p_flux = flux_posterior(xi, sd, pi, mu_f, t0, t1)
    p_state = joint_state_probs(ph["p_pheno"], p_flux)
    gauss_away = 1.0 - norm.cdf(0.0, loc=xi, scale=np.clip(sd, 1e-4, 10.0))
    return {
        "params": p,
        "xi_mean": xi,
        "p_away": p_flux[:, 1],
        "p_toward": p_flux[:, 2],
        "p_none": p_flux[:, 0],
        "p_flux": p_flux,
        "p_pheno": ph["p_pheno"],
        "p_state": p_state,
        "state": state_from_probs(p_state),
        "pheno": phenotype_from_probs(ph["p_pheno"]),
        "h": ph["h"],
        "h_als": ph["h_als"],
        "theta": ph["theta"],
        "beta": ph["beta"],
        "gauss_p_away": gauss_away,
        "xi_sd": sd,
        "loss": hist,
        "lambda": lam,
        "kappa": np.exp(np.clip(p["log_kappa"], -8, 8)),
        "rho_control": np.exp(np.clip(p["log_rho0"], -4, 4)),
        "c_s": p["c_s"],
        "c_g2m": p["c_g2m"],
        "sigma_v": sigma,
        "pi_flux": pi,
        "mu_flux": mu_f,
    }


def phenotype_calls(
    theta: np.ndarray,
    exposed: np.ndarray | None = None,
    *,
    force_control_reverted: bool = True,
) -> np.ndarray:
    """θ-gate labels. Forced control→reverted is a design constraint, not a measurement.

    Pass force_control_reverted=False for the empirical diagnostic (how many
    never-hypoxic cells would fail the persist/partial gates).
    """
    out = np.array(["partial"] * theta.size, dtype=object)
    out[theta <= 0.3] = "persistent"
    out[theta >= 0.7] = "reverted"
    if force_control_reverted and exposed is not None:
        out[np.asarray(exposed) < 0.5] = "reverted"
    return out


def state_calls(
    theta: np.ndarray, p_away: np.ndarray, p_toward: np.ndarray, exposed: np.ndarray
) -> np.ndarray:
    p_none = np.clip(1.0 - p_away - p_toward, 0.0, 1.0)
    p_flux = np.column_stack([p_none, p_away, p_toward])
    h_proxy = -np.log(np.clip(theta, 1e-6, 1 - 1e-6) / np.clip(1.0 - theta, 1e-6, None))
    p_pheno = phenotype_probs(h_proxy, exposed)
    return state_from_probs(joint_state_probs(p_pheno, p_flux))


def direction_calls(theta, p_away, p_toward, exposed):
    return state_calls(theta, p_away, p_toward, exposed)
