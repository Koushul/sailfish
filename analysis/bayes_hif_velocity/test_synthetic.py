#!/usr/bin/env python3
"""Fit the Bayesian lag model on a synthetic HIF-α target panel."""
from __future__ import annotations

import sys
from pathlib import Path

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
    print(f"AUROC transitioning_out vs partial: {m['auc_transition_vs_partial']:.3f}")
    print(
        f"theta persist={m['theta_persist_recall']:.3f} "
        f"partial={m['theta_partial_recall']:.3f} reverted={m['theta_reverted_recall']:.3f}"
    )
    print(
        f"recall partial={m['recall_partial']:.3f} t_out={m['recall_transitioning_out']:.3f} "
        f"t_in={m['recall_transitioning_in']:.3f}"
    )
    print(
        f"recall persist_exiting={m['recall_persistent_exiting']:.3f} "
        f"reverted_entering={m['recall_reverted_entering']:.3f}"
    )
    print(
        f"E[P] partial={m['prob_partial']:.3f} t_out={m['prob_t_out']:.3f} t_in={m['prob_t_in']:.3f}"
    )
    print(
        f"mean p_away  t_out={m['mean_p_away_reverting']:.3f}  "
        f"partial={m['mean_p_away_partial']:.3f}  control={m['mean_p_away_control']:.3f}"
    )
    print(f"mean p_toward t_in={m['mean_p_toward_inducing']:.3f}")
    print(f"mean |ξ| partial={m['mean_abs_xi_partial']:.3f} moving={m['mean_abs_xi_moving']:.3f}")
    print(f"control false flux rate: {m['fp_control']:.3f}")
    print(f"control persist (θ-gate): {m['ctrl_persist_theta']:.3f}")
    print(f"final loss {m['loss_final']:.1f}  start {m['loss_start']:.1f}")

    ok = True
    checks = [
        (m["spearman_xi"] > 0.32 and m["spearman_xi_p"] < 1e-6, "did not recover true lag"),
        (m["mean_p_away_reverting"] > m["mean_p_away_partial"] + 0.05, "t_out not higher p_away than partial"),
        (m["auc_transition_vs_partial"] > 0.7, "cannot separate transitioning_out from partial"),
        (m["auc_inducing"] > 0.65, "inducing AUROC too low"),
        (m["theta_partial_recall"] > 0.28, "GMM phenotype does not recover partial"),
        (m["recall_partial"] > 0.28, "partially reverted not recovered"),
        (m["prob_partial"] > 0.25, "posterior mass on partial too low"),
        (m["mean_abs_xi_moving"] > m["mean_abs_xi_partial"], "moving cells not larger |ξ| than partial"),
        (abs(m["spearman_cycle"]) < 0.35, "lag still tracks cell cycle"),
        (m["fp_control"] < 0.12, "control false flux rate too high"),
        (m["ctrl_persist_theta_forced"] == 0.0, "forced control labels are not all reverted"),
        (m["ctrl_persist_theta"] < 0.25, "empirical control persist tail too large (θ-gate overlaps never-hypoxic cells)"),
        (m["loss_final"] < m["loss_start"], "loss did not decrease"),
        (m["theta_persist_recall"] > 0.45 and m["theta_reverted_recall"] > 0.45, "persist vs reverted failed"),
        (m["auc_reverting"] > 0.68, "reverting AUROC too low"),
    ]
    for cond, msg in checks:
        if not cond:
            print(f"FAIL: {msg}")
            ok = False
    if ok:
        print("PASS")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
