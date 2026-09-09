from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np

from .io_utils import attach_adt, find_quants_h5ad, load_counts, prepare_adt, prepare_gex, write_mtx


def merge_gex_adt(
    *,
    gex_h5ad: Path | None,
    gex_alevin: Path | None,
    adt_h5ad: Path | None,
    adt_alevin: Path | None,
    feature_ref: Path | None,
    outdir: Path,
    sample: str,
    min_gex_umi: int,
    chemistry: str,
) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    gex = prepare_gex(load_counts(gex_h5ad, gex_alevin, "gex"))
    adt = None
    if adt_h5ad is not None or adt_alevin is not None:
        try:
            adt = prepare_adt(load_counts(adt_h5ad, adt_alevin, "adt"), feature_ref)
        except FileNotFoundError:
            adt = None

    keep = gex.obs["gex_counts"].values >= min_gex_umi
    gex_c = attach_adt(gex[keep].copy(), adt)
    gex_c.uns.update(
        {
            "sample": sample,
            "chemistry": chemistry,
            "mode": "quant",
            "cell_filter": f"GEX UMI>={min_gex_umi}; ADT left-join" if adt is not None else f"GEX UMI>={min_gex_umi}",
        }
    )
    out = outdir / f"{sample}_gex_adt.h5ad" if adt is not None else outdir / f"{sample}_gex.h5ad"
    gex_c.write_h5ad(out, compression="gzip")
    write_mtx(outdir / "filtered_feature_bc_matrix", gex_c)
    write_mtx(outdir / "raw_feature_bc_matrix", attach_adt(gex.copy(), adt))
    summary = {
        "sample": sample,
        "mode": "quant",
        "n_raw": int(gex.n_obs),
        "n_filtered": int(gex_c.n_obs),
        "n_vars": int(gex_c.n_vars),
        "median_gex_umi": float(np.median(gex_c.obs["gex_counts"])),
        "n_has_adt": int(gex_c.obs["has_adt"].sum()) if "has_adt" in gex_c.obs else 0,
        "min_gex_umi": min_gex_umi,
        "output": str(out),
    }
    (outdir / "quant_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return out
