#!/usr/bin/env python3
"""Sweep synthetic HIF-α panels: dropout, cycle, decoys, panel size, scrambled U."""
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
    ("core_default", SimConfig()),
    ("core_n8", SimConfig(n_core=8)),
    ("core_n16", SimConfig(n_core=16)),
    ("high_dropout", SimConfig(dropout=0.35)),
    ("strong_cycle", SimConfig(cycle_on_hif=0.45, cycle_on_counts=2.0)),
    ("no_capture_shift", SimConfig(capture_shift=0.0)),
    ("core_plus_decoy10", SimConfig(n_decoy=10)),
    ("weak_lag", SimConfig(lag_scale=0.35)),
    ("scrambled_U", SimConfig(scramble_unspliced=True)),
]


def run_one(name: str, cfg: SimConfig, seed: int) -> dict:
    cfg = SimConfig(**{**cfg.__dict__, "seed": seed})
    data, truth = simulate(cfg)
    est = fit(data, n_steps=400, lr=0.04, seed=seed + 17)
    m = evaluate_fit(data, truth, est)
    m["scenario"] = name
    m["seed"] = float(seed)
    return m


def main() -> int:
    rows: list[dict] = []
    seeds = (7, 11, 13)
    keys = None
    for name, cfg in SCENARIOS:
        for seed in seeds:
            m = run_one(name, cfg, seed)
            rows.append(m)
            print(
                f"{name:20s} seed={seed}  ρξ={m['spearman_xi']:.3f}  "
                f"cycle={m['spearman_cycle']:.3f}  AUCrev={m['auc_reverting']:.3f}  "
                f"AUCind={m['auc_inducing']:.3f}  θpersist={m['theta_persist_recall']:.3f}  "
                f"θrev={m['theta_reverted_recall']:.3f}  fp0={m['fp_control']:.3f}"
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

    by = {}
    for row in rows:
        by.setdefault(row["scenario"], []).append(row)

    def mean(name: str, key: str) -> float:
        xs = [r[key] for r in by[name]]
        return float(np_mean(xs))

    def np_mean(xs):
        import numpy as np

        return float(np.nanmean(np.asarray(xs, dtype=float)))

    ok = True
    if mean("core_default", "spearman_xi") < 0.35:
        print("FAIL: core panel lag recovery")
        ok = False
    if mean("core_default", "auc_reverting") < 0.65:
        print("FAIL: core reverting AUROC")
        ok = False
    if mean("scrambled_U", "spearman_xi") > mean("core_default", "spearman_xi") - 0.15:
        print("FAIL: scrambled unspliced should destroy lag recovery")
        ok = False
    if mean("core_plus_decoy10", "spearman_cycle") > 0.5:
        print("FAIL: decoy cycle genes leaked into ξ")
        ok = False
    if ok:
        print("BENCHMARK_PASS")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
