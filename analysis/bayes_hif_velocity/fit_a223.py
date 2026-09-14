#!/usr/bin/env python3
"""Fit E27 and E29 separately, then write a combined note. Does not write h5ads."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable
OUT = HERE / "results"


def main() -> None:
    for cfg in (HERE / "datasets" / "a223_e27.json", HERE / "datasets" / "a223_e29.json"):
        cmd = [PY, str(HERE / "fit_dataset.py"), "--config", str(cfg), "--min-control", "3"]
        print("+", " ".join(cmd), flush=True)
        subprocess.check_call(cmd, cwd=str(HERE.parent.parent))
    lines = [
        "# A223 E27 / E29 HIF-α lag",
        "",
        "`hypoxia_plus` is Image-iT LIVE Green ROS (DCF), not GFP and not HIF-α protein. Fits are **per chemistry** (E27 3′, E29 5′). Lineage comes from the tumor and neutrophil QC tables, not from every droplet in the OCM object. DN is the control arm; hypoxia_plus / DP / lactate_plus are exposed.",
        "",
        "DN n is below the recommended 40 in every lineage. Control MAD for \(h\) and the control mean of \(\\xi\) are poorly identified. Forced DN persist is 0% by design. Read empirical DN persist and per-gate exposed fractions. Do not transfer the E14 control location. Do not write into A223 h5ads.",
        "",
        "- E27: `a223_e27.md`",
        "- E29: `a223_e29.md`",
        "- Gate counts: `a223_gate_counts.tsv`",
        "",
    ]
    (OUT / "a223.md").write_text("\n".join(lines))
    print("wrote", OUT / "a223.md")


if __name__ == "__main__":
    main()
