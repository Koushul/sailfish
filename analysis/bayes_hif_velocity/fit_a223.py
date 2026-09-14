#!/usr/bin/env python3
"""Fit E27 and E29 separately, then write a combined note. Does not write h5ads."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
PY = sys.executable
OUT = HERE / "results"


def _stitch() -> None:
    rows = []
    for lane, summ_name, ev_name in (
        ("E27", "a223_e27_summary.tsv", "a223_e27_eval.tsv"),
        ("E29", "a223_e29_summary.tsv", "a223_e29_eval.tsv"),
    ):
        summ = pd.read_csv(OUT / summ_name, sep="\t")
        ev = pd.read_csv(OUT / ev_name, sep="\t")
        nctrl = {r.lineage: int(r.n_control) for r in ev.itertuples()}
        for r in summ.itertuples():
            rows.append(
                {
                    "lane": lane,
                    "lineage": r.lineage,
                    "ocm_gate": r.sample,
                    "n": int(r.n_cells),
                    "n_control_in_fit": nctrl.get(r.lineage),
                    "frac_persistent_forced": float(r.frac_pheno_persistent),
                    "frac_persistent_empirical": float(r.frac_pheno_empirical_persistent),
                    "mean_h": float(r.mean_h),
                    "mean_p_away": float(r.mean_p_away),
                }
            )
    tab = pd.DataFrame(rows)
    try:
        body = tab.to_markdown(index=False)
    except Exception:
        body = "```\n" + tab.to_string(index=False) + "\n```"
    lines = [
        "# A223 E27 / E29 HIF-α lag",
        "",
        "`hypoxia_plus` is Image-iT LIVE Green ROS (DCF), not GFP and not HIF-α protein. Fits are **per chemistry** (E27 3′, E29 5′). Lineage comes from the tumor and neutrophil QC tables, not from every droplet in the OCM object (E27 has 731 DN droplets; only 13 are QC tumors). DN is the control arm; hypoxia_plus / DP / lactate_plus are exposed. Do not transfer the E14 control location. Do not write into A223 h5ads.",
        "",
        "DN n is below 40 in every lineage. Control MAD for \(h\) and the control mean of \(\\xi\) are poorly identified. Forced DN persist is 0% by design. Neutrophil persist on E29 (~93% of exposed) is what a 3-cell control MAD does — not a HIF call. E27 vs E29 tumor persist (39% vs 7%) is not replicated across chemistry; UMI ratio is 0.52 vs 2.75. Lag \(p^{\\mathrm{away}}\) is ~0.01–0.02 in tumors for DN and DCF+ alike.",
        "",
        "## Per-gate fractions",
        "",
        body,
        "",
        "Per-library write-ups: `a223_e27.md`, `a223_e29.md`. Gate inventory: `a223_gate_counts.tsv`.",
        "",
    ]
    (OUT / "a223.md").write_text("\n".join(lines))
    tab.to_csv(OUT / "a223_fit_gates.tsv", sep="\t", index=False)
    print("wrote", OUT / "a223.md")


def main() -> None:
    for cfg in (HERE / "datasets" / "a223_e27.json", HERE / "datasets" / "a223_e29.json"):
        cmd = [PY, str(HERE / "fit_dataset.py"), "--config", str(cfg), "--min-control", "2"]
        print("+", " ".join(cmd), flush=True)
        subprocess.check_call(cmd, cwd=str(HERE.parent.parent))
    _stitch()


if __name__ == "__main__":
    main()
