#!/usr/bin/env python3
"""Shared metrics for synthetic HIF-α lag recovery."""
from __future__ import annotations

import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

from model import FLUX_STATES, STATE_NAMES


def roc_auc(y_true: np.ndarray, scores: np.ndarray) -> float:
    y = np.asarray(y_true, dtype=bool)
    s = np.asarray(scores, dtype=np.float64)
    n_pos = int(y.sum())
    n_neg = int((~y).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    u = mannwhitneyu(s[y], s[~y], alternative="greater").statistic
    return float(u / (n_pos * n_neg))


def _rate(pred: np.ndarray, mask: np.ndarray, label: str) -> float:
    if not np.any(mask):
        return float("nan")
    return float(np.mean(pred[mask] == label))


def _mean_prob(p_state: np.ndarray, mask: np.ndarray, label: str) -> float:
    if not np.any(mask) or p_state is None:
        return float("nan")
    j = STATE_NAMES.index(label)
    return float(p_state[mask, j].mean())


def evaluate_fit(data, truth: dict, est: dict) -> dict[str, float]:
    xi_hat = est["xi_mean"]
    xi = truth["xi"]
    rho, pval = spearmanr(xi_hat, xi)
    cycle_rho, _ = spearmanr(xi_hat, data.cycle_s)
    ctrl = data.exposed < 0.5
    kind = truth["kind"]
    persist = kind == "persist"
    partial = kind == "partial"
    reverted = kind == "reverted"
    t_out = kind == "transitioning_out"
    t_in = kind == "transitioning_in"
    p_exit = kind == "persistent_exiting"
    r_enter = kind == "reverted_entering"
    moving = t_out | t_in | p_exit | r_enter
    steady_exp = (~ctrl) & (kind == "partial")
    calls = est["state"]
    pheno = est["pheno"]
    p_state = est["p_state"]
    mid = partial | t_out

    return {
        "n_cells": float(data.U.shape[0]),
        "n_genes": float(data.U.shape[1]),
        "n_core": float(truth["n_core"]),
        "n_decoy": float(truth["n_decoy"]),
        "spearman_xi": float(rho),
        "spearman_xi_p": float(pval),
        "spearman_cycle": float(cycle_rho),
        "auc_reverting": roc_auc(t_out | p_exit, est["p_away"]),
        "auc_inducing": roc_auc(t_in | r_enter, est["p_toward"]),
        "auc_transition_vs_partial": roc_auc(t_out[mid], est["p_away"][mid]) if np.any(mid) else float("nan"),
        "auc_partial_theta": roc_auc(partial, est["p_pheno"][:, 1]),
        "mean_p_away_reverting": float(est["p_away"][t_out].mean()) if t_out.any() else float("nan"),
        "mean_p_away_partial": float(est["p_away"][partial].mean()) if partial.any() else float("nan"),
        "mean_p_away_control": float(est["p_away"][ctrl].mean()),
        "mean_p_toward_inducing": float(est["p_toward"][t_in].mean()) if t_in.any() else float("nan"),
        "fp_control": float(np.mean(np.isin(calls[ctrl], FLUX_STATES))),
        "recall_partial": _rate(calls, partial, "partial"),
        "recall_transitioning_out": _rate(calls, t_out, "transitioning_out"),
        "recall_transitioning_in": _rate(calls, t_in, "transitioning_in"),
        "recall_persistent_exiting": _rate(calls, p_exit, "persistent_exiting"),
        "recall_reverted_entering": _rate(calls, r_enter, "reverted_entering"),
        "recall_reverting": _rate(calls, t_out, "transitioning_out"),
        "recall_inducing": _rate(calls, t_in, "transitioning_in"),
        "prob_partial": _mean_prob(p_state, partial, "partial"),
        "prob_t_out": _mean_prob(p_state, t_out, "transitioning_out"),
        "prob_t_in": _mean_prob(p_state, t_in, "transitioning_in"),
        "theta_persist_recall": _rate(pheno, persist, "persistent"),
        "theta_partial_recall": _rate(pheno, partial, "partial"),
        "theta_reverted_recall": _rate(pheno, reverted, "reverted"),
        "mean_abs_xi_control": float(np.mean(np.abs(xi_hat[ctrl]))),
        "mean_abs_xi_partial": float(np.mean(np.abs(xi_hat[steady_exp]))) if steady_exp.any() else float("nan"),
        "mean_abs_xi_moving": float(np.mean(np.abs(xi_hat[moving]))) if moving.any() else float("nan"),
        "loss_start": float(est["loss"][0]),
        "loss_final": float(est["loss"][-1]),
    }
