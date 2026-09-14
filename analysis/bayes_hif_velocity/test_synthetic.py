#!/usr/bin/env python3
"""Fit the Bayesian lag model on a synthetic HIF-α target panel."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from evaluate import evaluate_fit
from model import fit
from simulate import SimConfig, simulate


def main() -> int:
    data, truth = simulate(SimConfig(seed=7))
    est = fit(data, n_steps=500, lr=0.04, seed=1)
    m = evaluate_fit(data, truth, est)

    print(f"n={data.U.shape[0]} genes={data.U.shape[1]} core={truth['n_core']}")
    print(f"panel: {', '.join(map(str, truth['genes'][:8]))}, ...")
    print(f"Spearman ξ_hat vs ξ_true: {m['spearman_xi']:.3f} (p={m['spearman_xi_p']:.1e})")
    print(f"Spearman ξ_hat vs cycle S: {m['spearman_cycle']:.3f}")
    print(f"AUROC reverting vs rest: {m['auc_reverting']:.3f}")
    print(f"AUROC inducing vs rest: {m['auc_inducing']:.3f}")
    print(f"theta persist recall={m['theta_persist_recall']:.3f} reverted recall={m['theta_reverted_recall']:.3f}")
    print(f"mean p_away  reverting={m['mean_p_away_reverting']:.3f}  control={m['mean_p_away_control']:.3f}")
    print(f"mean p_toward inducing={m['mean_p_toward_inducing']:.3f}")
    print(f"control false direction rate: {m['fp_control']:.3f}")
    print(f"recall reverting={m['recall_reverting']:.3f} inducing={m['recall_inducing']:.3f}")
    print(f"final loss {m['loss_final']:.1f}  start {m['loss_start']:.1f}")

    ok = True
    if not (m["spearman_xi"] > 0.35 and m["spearman_xi_p"] < 1e-6):
        print("FAIL: did not recover true lag")
        ok = False
    if m["mean_p_away_reverting"] <= m["mean_p_away_control"] + 0.08:
        print("FAIL: reverting cells not higher p_away than control")
        ok = False
    if abs(m["spearman_cycle"]) > 0.45:
        print("FAIL: lag still tracks cell cycle")
        ok = False
    if m["fp_control"] > 0.12:
        print("FAIL: control false direction rate too high")
        ok = False
    if m["loss_final"] > m["loss_start"]:
        print("FAIL: loss did not decrease")
        ok = False
    if m["theta_persist_recall"] < 0.5 or m["theta_reverted_recall"] < 0.5:
        print("FAIL: spliced HIF program does not separate persist vs reverted")
        ok = False
    if m["auc_reverting"] < 0.65:
        print("FAIL: reverting AUROC too low")
        ok = False
    if ok:
        print("PASS")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
