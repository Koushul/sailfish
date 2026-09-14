#!/usr/bin/env python3
"""Shared metrics for synthetic HIF-α lag recovery."""
from __future__ import annotations

import numpy as np
from scipy.stats import mannwhitneyu, spearmanr

from model import direction_calls, phenotype_calls


def roc_auc(y_true: np.ndarray, scores: np.ndarray) -> float:
    y = np.asarray(y_true, dtype=bool)
    s = np.asarray(scores, dtype=np.float64)
    n_pos = int(y.sum())
    n_neg = int((~y).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    u = mannwhitneyu(s[y], s[~y], alternative="greater").statistic
    return float(u / (n_pos * n_neg))


def evaluate_fit(data, truth: dict, est: dict) -> dict[str, float]:
    xi_hat = est["xi_mean"]
    xi = truth["xi"]
    rho, pval = spearmanr(xi_hat, xi)
    cycle_rho, _ = spearmanr(xi_hat, data.cycle_s)
    ctrl = data.exposed < 0.5
    exp = ~ctrl
    kind = truth["kind"]
    rev = kind == "reverting"
    ind = kind == "inducing"
    persist = kind == "persist"
    reverted = kind == "reverted"
    rest_exp = exp & ~rev & ~ind

    calls = direction_calls(data.theta, est["p_away"], est["p_toward"], data.exposed)
    pheno = phenotype_calls(data.theta)

    return {
        "n_cells": float(data.U.shape[0]),
        "n_genes": float(data.U.shape[1]),
        "n_core": float(truth["n_core"]),
        "n_decoy": float(truth["n_decoy"]),
        "spearman_xi": float(rho),
        "spearman_xi_p": float(pval),
        "spearman_cycle": float(cycle_rho),
        "auc_reverting": roc_auc(rev, est["p_away"]),
        "auc_inducing": roc_auc(ind, est["p_toward"]),
        "mean_p_away_reverting": float(est["p_away"][rev].mean()) if rev.any() else float("nan"),
        "mean_p_away_control": float(est["p_away"][ctrl].mean()),
        "mean_p_toward_inducing": float(est["p_toward"][ind].mean()) if ind.any() else float("nan"),
        "fp_control": float(np.mean(np.isin(calls[ctrl], ["reverting", "inducing"]))),
        "recall_reverting": float(np.mean(calls[rev] == "reverting")) if rev.any() else float("nan"),
        "recall_inducing": float(np.mean(calls[ind] == "inducing")) if ind.any() else float("nan"),
        "theta_persist_recall": float(np.mean(pheno[persist] == "persistent")) if persist.any() else float("nan"),
        "theta_reverted_recall": float(np.mean(pheno[reverted] == "reverted")) if reverted.any() else float("nan"),
        "mean_abs_xi_control": float(np.mean(np.abs(xi_hat[ctrl]))),
        "mean_abs_xi_steady_exposed": float(np.mean(np.abs(xi_hat[rest_exp]))) if rest_exp.any() else float("nan"),
        "loss_start": float(est["loss"][0]),
        "loss_final": float(est["loss"][-1]),
    }
