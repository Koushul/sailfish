#!/usr/bin/env python3
"""Fit the HIF-α lag model on A223 E27 tumor cells only (single cohort).

Uses the E27 OCM object and tumor QC barcodes. Does not use E14S, E15S, E29,
neutrophils, or DN as a control arm. Image-iT gates are recorded after the fit.
Does not write into the h5ad.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from fit_dataset import overlay_qc_tables
from load_h5ad import load_counts
from single_cohort import LINEAGE, MIN_SPLICED, add_args, run_cohort

H5AD = "/ix1/ylee/kor11/A223/E27/ocm/E27_gex_adt_ocm.h5ad"
TUMOR_QC = "/ix1/ylee/kor11/A223/tumor_kinetics/tumor_qc.csv"


def main() -> None:
    p = argparse.ArgumentParser()
    add_args(p)
    p.set_defaults(h5ad=H5AD)
    args = p.parse_args()
    placed = load_counts(
        args.h5ad,
        sample_col="sample_id",
        lineage_col=None,
        gene_symbol_col="gene_symbol",
        historical_state_col=None,
        historical_theta_col=None,
        require_all_panel=False,
    )
    overlay_qc_tables(placed, [{"path": TUMOR_QC, "lane": "E27", "lineage": "Tumor"}])
    keep = (placed.cell_group == LINEAGE) & (placed.L >= MIN_SPLICED)
    intro = (
        "A223 E27 (3′ OCM GEX), lineage `Tumor`, spliced UMI ≥ 5000. "
        "Tumor barcodes come from the tumor QC table, not from every droplet. "
        "No E14S, E15S, or E29 cells. No neutrophil pool. "
        "Image-iT gates (`hypoxia_plus` is the same GFP-channel probe as E15S; `DN` is Image-iT−) "
        "are **not** used as control vs exposed. "
        r"\(h\) is this cohort's own median/MAD among E27 tumors. "
        "Persist means ≥1.5 MAD above a typical E27 tumor in this fit. "
        "Do not write into the h5ad."
    )
    run_cohort(
        placed,
        keep,
        slug="e27",
        cohort_name="E27",
        args=args,
        extra_cell_cols={"ocm_gate": placed.sample},
        intro=intro,
    )


if __name__ == "__main__":
    main()
