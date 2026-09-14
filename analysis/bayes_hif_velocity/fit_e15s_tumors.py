#!/usr/bin/env python3
"""Fit the HIF-α lag model on E15S tumor cells only (single cohort).

Does not use E14S cells or any parameters estimated on E14S.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from single_cohort import add_args, run


def main() -> None:
    p = argparse.ArgumentParser()
    add_args(p)
    run("E15S", "e15s", p.parse_args())


if __name__ == "__main__":
    main()
