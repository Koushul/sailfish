#!/usr/bin/env python3
"""Build the 2-oligo E14S/E15S localization site (data.js + cells.js)."""
from __future__ import annotations

import argparse
import base64
import json
import re
import shutil
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from cell_placement import (  # noqa: E402
    BETA,
    GRID,
    MAX_ENTROPY_BITS,
    axis_prior,
    compact_oligo_sequence,
    extract_core,
    peak_norm_u8,
    place_from_axis_counts,
    write_data_js,
)

CHIP_LOC_RE = re.compile(r"^(ROW|COLUMN)\s+(\d+)$", re.I)


def layout_barcodes_from_csv(path: Path) -> list[dict]:
    df = pd.read_csv(path)
    df = df.dropna(subset=["Name", "Chip Location (Row/Column)"])
    barcodes = []
    for _, row in df.iterrows():
        loc = str(row["Chip Location (Row/Column)"]).strip()
        m = CHIP_LOC_RE.match(loc)
        if not m:
            raise SystemExit(f"Unparseable chip location {loc!r} for {row['Name']}")
        axis = "row" if m.group(1).upper() == "ROW" else "column"
        plate = int(str(row["96 Source Plate #"]).replace("Plate", "").strip())
        barcodes.append(
            {
                "name": str(row["Name"]),
                "plate": plate,
                "axis": axis,
                "index": int(m.group(2)),
                "well": str(row["Well Position in 384 Source Plate"]),
                "sequence": compact_oligo_sequence(row["Sequence"]),
            }
        )
    return barcodes


def load_two_oligo_layout(path: Path) -> tuple[list[dict], list[dict]]:
    rows: dict[int, dict] = {}
    cols: dict[int, dict] = {}
    for b in layout_barcodes_from_csv(path):
        item = {**b, "core": extract_core(b["sequence"])}
        if b["plate"] == 1 and b["axis"] == "row":
            rows[b["index"]] = item
        elif b["plate"] == 2 and b["axis"] == "column":
            cols[b["index"]] = item
    missing_r = [i for i in range(1, GRID + 1) if i not in rows]
    missing_c = [i for i in range(1, GRID + 1) if i not in cols]
    if missing_r or missing_c:
        raise SystemExit(f"layout missing row indices {missing_r} col indices {missing_c}")
    return [rows[i] for i in range(1, GRID + 1)], [cols[i] for i in range(1, GRID + 1)]


def adt_axis_counts(adt: np.ndarray, adt_names: list[str], barcodes: list[dict], feature_ref: Path) -> np.ndarray:
    ref = pd.read_csv(feature_ref)
    name_to_seq = dict(zip(ref["name"].astype(str), ref["sequence"].astype(str)))
    seq_to_col = {}
    for j, n in enumerate(adt_names):
        seq = name_to_seq.get(str(n))
        if seq:
            seq_to_col[seq] = j
    cols = []
    for bc in barcodes:
        if bc["core"] not in seq_to_col:
            raise SystemExit(f"No ADT feature for {bc['name']} seq={bc['core']}")
        cols.append(seq_to_col[bc["core"]])
    return np.asarray(adt, dtype=np.float64)[:, cols]


def b64_arr(a: np.ndarray) -> str:
    return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode("ascii")


def f32_list(a: np.ndarray) -> list[float]:
    return [float(x) for x in np.asarray(a, dtype=np.float64).ravel()]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--h5ad",
        default="/ix1/ylee/kor11/MC38/E14SE15S/E14SE15S_gex_adt_placed.h5ad",
    )
    ap.add_argument("--layout", default=str(HERE / "layouts" / "e14se15s_layout.csv"))
    ap.add_argument(
        "--feature-ref",
        default=str(HERE.parent / "refs" / "new_feature_ref_quant.csv"),
    )
    ap.add_argument("--site", default=str(HERE / "sites" / "e14se15s_2oligo"))
    args = ap.parse_args()

    site = Path(args.site)
    cells_dir = site / "datasets" / "E14SE15S"
    cells_dir.mkdir(parents=True, exist_ok=True)

    barcodes = layout_barcodes_from_csv(Path(args.layout))
    write_data_js(barcodes, site / "data.js")

    row_bc, col_bc = load_two_oligo_layout(Path(args.layout))
    adata = ad.read_h5ad(args.h5ad)
    adt = adata.obsm["ADT"]
    adt = adt.toarray() if sparse.issparse(adt) else np.asarray(adt, dtype=np.float64)
    adt_names = list(adata.uns["ADT_var"]["feature_name"])
    row_counts = adt_axis_counts(adt, adt_names, row_bc, Path(args.feature_ref))
    col_counts = adt_axis_counts(adt, adt_names, col_bc, Path(args.feature_ref))
    placed = place_from_axis_counts(row_counts, col_counts)

    obs_map_r = adata.obs["map_row"].to_numpy(dtype=int)
    obs_map_c = adata.obs["map_col"].to_numpy(dtype=int)
    n_mismatch = int(np.sum((placed["map_row"] != obs_map_r) | (placed["map_col"] != obs_map_c)))
    if n_mismatch:
        print(f"warning: {n_mismatch} cells differ from h5ad MAP wells (using recomputed)")

    sample = adata.obs["sample"].astype(str).tolist()
    payload = {
        "layout": "two_oligo",
        "grid_size": GRID,
        "n_cells": int(adata.n_obs),
        "max_entropy_bits": MAX_ENTROPY_BITS,
        "beta": BETA,
        "obs_names": list(map(str, adata.obs_names)),
        "sample_id": sample,
        "map_row": placed["map_row"].astype(int).tolist(),
        "map_col": placed["map_col"].astype(int).tolist(),
        "confidence": f32_list(placed["confidence"]),
        "row_entropy": f32_list(placed["row_entropy"]),
        "col_entropy": f32_list(placed["col_entropy"]),
        "total_entropy": f32_list(placed["total_entropy"]),
        "spatial_entropy": f32_list(placed["spatial_entropy"]),
        "total_counts": placed["layout_umi"].astype(int).tolist(),
        "row_post_b64": b64_arr(peak_norm_u8(placed["row_post"])),
        "col_post_b64": b64_arr(peak_norm_u8(placed["col_post"])),
        "row_ll_b64": b64_arr(np.ascontiguousarray(placed["row_ll"], dtype=np.float32)),
        "col_ll_b64": b64_arr(np.ascontiguousarray(placed["col_ll"], dtype=np.float32)),
        "row_prior": f32_list(axis_prior(row_counts)),
        "col_prior": f32_list(axis_prior(col_counts)),
    }
    cells_path = cells_dir / "cells.js"
    cells_path.write_text("window.CELL_DATA = " + json.dumps(payload, separators=(",", ":")) + ";\n")
    print(f"wrote {cells_path} n={adata.n_obs} bytes={cells_path.stat().st_size}")

    csv_src = Path(args.h5ad).with_suffix(".assignments.csv")
    csv_dst = site / "assignments.csv"
    if csv_src.exists():
        shutil.copy2(csv_src, csv_dst)
        print(f"copied {csv_src} -> {csv_dst}")
    else:
        pd.DataFrame(
            {
                "barcode": adata.obs_names,
                "sample_id": sample,
                "map_row": placed["map_row"],
                "map_col": placed["map_col"],
                "confidence": placed["confidence"],
                "row_entropy": placed["row_entropy"],
                "col_entropy": placed["col_entropy"],
                "total_entropy": placed["total_entropy"],
                "spatial_entropy": placed["spatial_entropy"],
                "layout_umi": placed["layout_umi"],
            }
        ).to_csv(csv_dst, index=False)
        print(f"wrote {csv_dst}")

    print(
        f"wells={pd.Series(list(zip(placed['map_row'], placed['map_col']))).nunique()} "
        f"median H={np.median(placed['total_entropy']):.3f} "
        f"samples={dict(pd.Series(sample).value_counts())}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
