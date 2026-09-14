#!/usr/bin/env python3
"""Sweep synthetic HIF-α panels under realistic confounders."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from evaluate import evaluate_fit
from model import fit
from simulate import SimConfig, simulate

SCENARIOS: list[tuple[str, SimConfig]] = [
    ("realistic", SimConfig()),
    ("core_n8", SimConfig(n_core=8, n_decoy=0, n_silent=0)),
    ("high_dropout", SimConfig(dropout=0.4)),
    ("strong_cycle", SimConfig(cycle_on_hif=0.45, cycle_on_spliced=0.55, cycle_on_counts=1.8)),
    ("weak_lag", SimConfig(lag_scale=0.4)),
    ("no_capture_shift", SimConfig(capture_shift=0.0)),
    ("scrambled_U", SimConfig(scramble_unspliced=True)),
]


def run_one(name: str, cfg: SimConfig, seed: int) -> dict:
    cfg = SimConfig(**{**cfg.__dict__, "seed": seed})
    data, truth = simulate(cfg)
    est = fit(data, n_steps=420, lr=0.04, seed=seed + 17)
    m = evaluate_fit(data, truth, est)
    m["scenario"] = name
    m["seed"] = float(seed)
    return m


def main() -> int:
    import numpy as np

    rows: list[dict] = []
    keys = None
    for name, cfg in SCENARIOS:
        for seed in (7, 11, 13):
            m = run_one(name, cfg, seed)
            rows.append(m)
            print(
                f"{name:18s} seed={seed}  ρξ={m['spearman_xi']:.3f}  "
                f"AUCrev={m['auc_reverting']:.3f} AUCind={m['auc_inducing']:.3f}  "
                f"t_vs_p={m['auc_transition_vs_partial']:.3f}  "
                f"part={m['recall_partial']:.3f} t_out={m['recall_transitioning_out']:.3f} "
                f"t_in={m['recall_transitioning_in']:.3f}  "
                f"θpart={m['theta_partial_recall']:.3f} fp0={m['fp_control']:.3f}"
            )
            keys = list(m.keys())

    out = HERE / "results"
    out.mkdir(exist_ok=True)
    path = out / "synthetic_benchmark.tsv"
    fieldnames = ["scenario", "seed"] + [k for k in keys if k not in ("scenario", "seed")]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow({k: row[k] for k in fieldnames})
    print(f"wrote {path}")

    by: dict[str, list] = {}
    for row in rows:
        by.setdefault(row["scenario"], []).append(row)

    def mean(name: str, key: str) -> float:
        return float(np.nanmean(np.asarray([r[key] for r in by[name]], dtype=float)))

    ok = True
    if mean("realistic", "spearman_xi") < 0.22:
        print("FAIL: realistic lag recovery")
        ok = False
    if mean("realistic", "auc_transition_vs_partial") < 0.85:
        print("FAIL: realistic t_out vs partial")
        ok = False
    if mean("realistic", "recall_partial") < 0.7:
        print("FAIL: realistic partial recall")
        ok = False
    if mean("realistic", "auc_inducing") < 0.65:
        print("FAIL: realistic inducing AUROC")
        ok = False
    if mean("realistic", "theta_partial_recall") < 0.85:
        print("FAIL: realistic partial phenotype")
        ok = False
    if mean("strong_cycle", "theta_partial_recall") < 0.7:
        print("FAIL: strong cycle swallowed partial phenotype")
        ok = False
    if mean("strong_cycle", "auc_transition_vs_partial") < 0.7:
        print("FAIL: strong cycle destroyed transition vs partial")
        ok = False
    if mean("scrambled_U", "auc_transition_vs_partial") > 0.7:
        print("FAIL: scrambled unspliced still separates transition")
        ok = False
    if mean("scrambled_U", "spearman_xi") > 0.15:
        print("FAIL: scrambled unspliced still recovers lag")
        ok = False
    if ok:
        print("BENCHMARK_PASS")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
