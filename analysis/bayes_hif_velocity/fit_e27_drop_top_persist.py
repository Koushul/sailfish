#!/usr/bin/env python3
"""E27 tumor control: drop the 600 highest-h cells from the previous fit and refit.

The second fit does not see those cells and does not use their parameters.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from fit_dataset import overlay_qc_tables
from load_h5ad import load_counts
from single_cohort import HERE, LINEAGE, MIN_SPLICED, add_args, run_cohort

H5AD = "/ix1/ylee/kor11/A223/E27/ocm/E27_gex_adt_ocm.h5ad"
TUMOR_QC = "/ix1/ylee/kor11/A223/tumor_kinetics/tumor_qc.csv"
PREV_CELLS = HERE / "results" / "e27_tumor_cells.tsv"
N_DROP = 600
SLUG = "e27_drop600"


def main() -> None:
    p = argparse.ArgumentParser()
    add_args(p)
    p.set_defaults(h5ad=H5AD)
    args = p.parse_args()
    prev = pd.read_csv(PREV_CELLS, sep="\t")
    if len(prev) <= N_DROP:
        raise SystemExit(f"previous E27 table has {len(prev)} cells; cannot drop {N_DROP}")
    removed = prev.nlargest(N_DROP, "h").copy()
    drop_ids = set(removed["cell"].astype(str))
    n_persist_removed = int((removed["pheno"] == "persistent").sum())
    h_min = float(removed["h"].min())
    h_max = float(removed["h"].max())

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
    ids = np.array(placed.cell_id, dtype=object).astype(str)
    keep = (placed.cell_group == LINEAGE) & (placed.L >= MIN_SPLICED) & ~np.isin(ids, list(drop_ids))
    n_keep = int(keep.sum())
    if n_keep < 20:
        raise SystemExit(f"too few cells after drop: {n_keep}")
    intro = (
        "Control: A223 E27 tumors after deleting the 600 highest-\(h\) cells from the previous "
        f"single-cohort fit (`e27_tumor_cells.tsv`). Removed \(h\) range {h_min:.2f}–{h_max:.2f}; "
        f"{n_persist_removed} of those 600 were θ-gate persistent in the first fit "
        f"(first-fit persist n=431). Remaining n={n_keep}. "
        "The model is fit from scratch on the remaining cells only: new gene detection, new ALS \(h\), "
        "new MAD, new lag. It is not told that any cells were removed. "
        "Image-iT gates still do not enter the fit. Do not write into the h5ad."
    )
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    removed.to_csv(out / f"{SLUG}_removed.tsv", sep="\t", index=False)
    run_cohort(
        placed,
        keep,
        slug=SLUG,
        cohort_name="E27 drop-top-600",
        args=args,
        extra_cell_cols={"ocm_gate": placed.sample},
        intro=intro,
    )

    new = pd.read_csv(out / f"{SLUG}_tumor_cells.tsv", sep="\t")
    lines = (out / f"{SLUG}_tumors.md").read_text().rstrip() + "\n"
    lines += "\n## Versus the untrimmed E27 fit\n\n"
    lines += (
        f"- Removed 600 cells with highest first-fit \(h\) (min removed \(h\)={h_min:.3f}). "
        f"Overlap of remaining barcodes with the removed list: "
        f"{int(np.isin(new['cell'].astype(str), list(drop_ids)).sum())}.\n"
        f"- First-fit persist fraction among all 5066: {float((prev['pheno']=='persistent').mean()):.3f}. "
        f"This refit's persist fraction: {float((new['pheno']=='persistent').mean()):.3f}.\n"
        "- If persist is only the upper tail of a unimodal factor, trimming the old tail and "
        "re-standardizing should grow a new tail of similar size.\n"
    )
    (out / f"{SLUG}_tumors.md").write_text(lines)


if __name__ == "__main__":
    main()
