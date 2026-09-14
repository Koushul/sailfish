#!/usr/bin/env python3
"""Fit the Bayesian lag model on synthetic counts and check recovery."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from model import direction_calls, fit
from simulate import simulate


def main() -> int:
    data, truth = simulate(seed=7)
    est = fit(data, n_steps=500, lr=0.04, seed=1)
    xi_hat = est["xi_mean"]
    xi = truth["xi"]
    rho, pval = spearmanr(xi_hat, xi)
    ctrl = data.exposed < 0.5
    exp = ~ctrl
    rev = truth["kind"] == "reverting"
    ind = truth["kind"] == "inducing"
    rest = exp & ~rev & ~ind

    calls = direction_calls(data.theta, est["p_away"], est["p_toward"], data.exposed)
    fp_ctrl = float(np.mean(np.isin(calls[ctrl], ["reverting", "inducing"])))
    rec_rev = float(np.mean(calls[rev] == "reverting")) if rev.any() else 0.0
    rec_ind = float(np.mean(calls[ind] == "inducing")) if ind.any() else 0.0
    cycle_rho, _ = spearmanr(xi_hat, data.cycle_s)
    mean_away_rev = float(est["p_away"][rev].mean()) if rev.any() else 0.0
    mean_away_ctrl = float(est["p_away"][ctrl].mean())
    mean_toward_ind = float(est["p_toward"][ind].mean()) if ind.any() else 0.0

    print(f"n={data.U.shape[0]} genes={data.U.shape[1]}")
    print(f"Spearman ξ_hat vs ξ_true: {rho:.3f} (p={pval:.1e})")
    print(f"Spearman ξ_hat vs cycle S: {cycle_rho:.3f}")
    print(f"mean p_away  reverting={mean_away_rev:.3f}  control={mean_away_ctrl:.3f}")
    print(f"mean p_toward inducing={mean_toward_ind:.3f}")
    print(f"control false direction rate: {fp_ctrl:.3f}")
    print(f"recall reverting={rec_rev:.3f} inducing={rec_ind:.3f}")
    print(f"final loss {est['loss'][-1]:.1f}  start {est['loss'][0]:.1f}")
    print(f"mean |ξ| control={np.mean(np.abs(xi_hat[ctrl])):.3f} rest_exposed={np.mean(np.abs(xi_hat[rest])):.3f}")

    ok = True
    if not (rho > 0.35 and pval < 1e-6):
        print("FAIL: did not recover true lag")
        ok = False
    if mean_away_rev <= mean_away_ctrl + 0.08:
        print("FAIL: reverting cells not higher p_away than control")
        ok = False
    if abs(cycle_rho) > 0.45:
        print("FAIL: lag still tracks cell cycle")
        ok = False
    if fp_ctrl > 0.12:
        print("FAIL: control false direction rate too high")
        ok = False
    if est["loss"][-1] > est["loss"][0]:
        print("FAIL: loss did not decrease")
        ok = False
    if ok:
        print("PASS")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
