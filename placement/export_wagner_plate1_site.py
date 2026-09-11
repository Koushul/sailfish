#!/usr/bin/env python3
"""Place MC38 cells and write a lucid-crystal-style site.

E15S GEX barcodes come from Palak Cell Ranger
(`/ix1/ylee/Palak/cellranger_apps/e15s/outs`), which is GEX-only.
Spatial-hash ADT is taken from the original E15S Cell Ranger raw matrix
(Antibody Capture), which covers every Palak barcode.

E14S ADT comes from the original E14S filtered matrix. Localization uses
both spatial-hash plates (the 48×48 chip needs Plate 1 rows + Plate 2 cols).
"""
from __future__ import annotations

import argparse
import base64
import gzip
import json
import shutil
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from cell_placement import (  # noqa: E402
    BETA,
    GRID,
    MAX_ENTROPY_BITS,
    axis_prior,
    extract_core,
    load_data_js,
    peak_norm_u8,
    place_from_counts,
)

DEFAULT_H5AD = Path("/ix1/ylee/shared/external/data/WagnerCollab/mc38_velocity.h5ad")
PALAK_E15S_BCS = Path(
    "/ix1/ylee/Palak/cellranger_apps/e15s/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
)
E15S_ADT_H5 = Path("/ix1/ylee/shared/MC38_Hypoxia_001/E15S/raw_feature_bc_matrix.h5")
E14S_ADT_H5 = Path("/ix1/ylee/shared/MC38_Hypoxia_001/E14S/filtered_feature_bc_matrix.h5")
DEFAULT_FEATURE_REF = HERE.parent / "refs" / "new_feature_ref_quant.csv"
DEFAULT_DATA_JS = HERE / "layouts" / "data.js"
DEFAULT_TEMPLATE = HERE / "site_template" / "index.html"
MIN_LAYOUT_UMI = 10


def b64_f32(arr: np.ndarray) -> str:
    return base64.b64encode(np.ascontiguousarray(arr, dtype=np.float32).tobytes()).decode("ascii")


def b64_u8(arr: np.ndarray) -> str:
    return base64.b64encode(np.ascontiguousarray(arr, dtype=np.uint8).tobytes()).decode("ascii")


def f32_list(x: np.ndarray, nd: int = 6) -> list[float]:
    return [round(float(v), nd) for v in np.asarray(x, dtype=np.float64)]


def _as_str(arr: np.ndarray) -> np.ndarray:
    out = []
    for x in arr:
        if isinstance(x, (bytes, np.bytes_)):
            out.append(x.decode())
        else:
            out.append(str(x))
    return np.array(out, dtype=object)


def adt_from_10x_h5(h5_path: Path, want: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Slice Antibody Capture columns for `want` barcodes from a 10x feature-barcode h5."""
    with h5py.File(h5_path, "r") as f:
        m = f["matrix"]
        barcodes = _as_str(m["barcodes"][:])
        lookup = {b: i for i, b in enumerate(barcodes)}
        ftype = _as_str(m["features"]["feature_type"][:])
        feat_names = _as_str(m["features"]["name"][:])
        adt_rows = np.where(ftype == "Antibody Capture")[0]
        if adt_rows.size == 0:
            raise SystemExit(f"No Antibody Capture features in {h5_path}")
        names = feat_names[adt_rows]
        row_to_k = {int(r): k for k, r in enumerate(adt_rows)}
        n_adt = int(adt_rows.size)
        n = len(want)
        out = np.zeros((n, n_adt), dtype=np.float64)
        found = np.zeros(n, dtype=bool)
        indptr = m["indptr"]
        indices = m["indices"]
        data = m["data"]
        for i, bc in enumerate(want):
            j = lookup.get(str(bc))
            if j is None:
                j = lookup.get(str(bc).split("-")[0] + "-1")
            if j is None:
                continue
            s, e = int(indptr[j]), int(indptr[j + 1])
            feat_idx = indices[s:e]
            vals = data[s:e]
            for fi, v in zip(feat_idx, vals):
                k = row_to_k.get(int(fi))
                if k is not None:
                    out[i, k] = float(v)
            found[i] = True
        print(f"ADT from {h5_path.name}: matched {int(found.sum())} / {n}")
        return out, names, found


def plate_blocks_from_adt(adt: np.ndarray, adt_names: np.ndarray, data_js: str, feature_ref: str):
    layout = load_data_js(data_js)
    ref = pd.read_csv(feature_ref)
    name_to_seq = dict(zip(ref["name"].astype(str), ref["sequence"].astype(str)))
    seq_to_col = {}
    for j, n in enumerate(adt_names):
        seq = name_to_seq.get(str(n))
        if seq:
            seq_to_col[seq] = j
    blocks = {}
    for axis_js, axis_key in (("row", "row"), ("column", "col")):
        for plate in (1, 2):
            cols = []
            for idx in range(1, GRID + 1):
                bc = layout["by"][(plate, axis_js, idx)]
                core = extract_core(bc["sequence"])
                if core not in seq_to_col:
                    raise SystemExit(f"No ADT feature for {bc['name']} seq={core}")
                cols.append(seq_to_col[core])
            blocks[f"{axis_key}_p{plate}"] = adt[:, cols].astype(np.float64)
    return blocks


def write_cells_js(path: Path, names: np.ndarray, sample: np.ndarray, out: dict, blocks: dict) -> None:
    payload = {
        "grid_size": GRID,
        "n_cells": int(len(names)),
        "max_entropy_bits": MAX_ENTROPY_BITS,
        "beta": BETA,
        "obs_names": names.tolist(),
        "sample_id": sample.tolist(),
        "map_row": out["map_row"].astype(int).tolist(),
        "map_col": out["map_col"].astype(int).tolist(),
        "confidence": f32_list(out["confidence"]),
        "row_entropy": f32_list(out["row_entropy"]),
        "col_entropy": f32_list(out["col_entropy"]),
        "total_entropy": f32_list(out["total_entropy"]),
        "spatial_entropy": f32_list(out["spatial_entropy"]),
        "total_counts": out["layout_umi"].astype(int).tolist(),
        "row_post_b64": b64_u8(peak_norm_u8(out["row_post"])),
        "col_post_b64": b64_u8(peak_norm_u8(out["col_post"])),
        "row_ll_p1_b64": b64_f32(out["row_ll_p1"]),
        "row_ll_p2_b64": b64_f32(out["row_ll_p2"]),
        "col_ll_p1_b64": b64_f32(out["col_ll_p1"]),
        "col_ll_p2_b64": b64_f32(out["col_ll_p2"]),
        "row_prior_p1": f32_list(axis_prior(blocks["row_p1"])),
        "row_prior_p2": f32_list(axis_prior(blocks["row_p2"])),
        "col_prior_p1": f32_list(axis_prior(blocks["col_p1"])),
        "col_prior_p2": f32_list(axis_prior(blocks["col_p2"])),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("window.CELL_DATA = " + json.dumps(payload, separators=(",", ":")) + ";\n")
    print(f"wrote {path} ({path.stat().st_size / 1e6:.2f} MB, n={len(names)})")


def patch_html(src: Path, n_cells: int) -> str:
    html = src.read_text()
    html = html.replace(
        """        <select id="dataset-sel">
          <option value="E28S" selected>E28S</option>
          <option value="E14S">E14S</option>
          <option value="E15S">E15S</option>
        </select>""",
        """        <select id="dataset-sel">
          <option value="MC38" selected>MC38 E14S+E15S</option>
        </select>""",
    )
    html = html.replace(
        """    const DATASETS = {
      E28S: { path: 'datasets/E28S/cells.js', note: 'edge/core hypoxia OCM · original published assignment' },
      E14S: { path: 'datasets/E14S/cells.js', note: 'entropy assignment · layout UMI≥10, GEX≥500' },
      E15S: { path: 'datasets/E15S/cells.js', note: 'entropy assignment · layout UMI≥10, GEX≥500' },
    };""",
        """    const DATASETS = {
      MC38: { path: 'datasets/MC38/cells.js', note: 'E15S GEX: Palak cellranger e15s/outs · ADT: original E15S raw · both plates · layout UMI≥10' },
    };""",
    )
    html = html.replace("let currentDataset = 'E28S';", "let currentDataset = 'MC38';")
    html = html.replace(
        "Hover a microwell to inspect barcodes; switch to Cells to position all cells and assess entropy",
        "MC38 · E15S from Palak cellranger + original ADT · both plates · "
        f"{n_cells:,} cells (layout UMI≥{MIN_LAYOUT_UMI}) · "
        "<a href='assignments.csv' style='color:#6bcf9e'>assignments.csv</a>",
    )
    html = html.replace(
        '<meta property="og:url" content="https://lucid-crystal-kmqy.here.now/" />',
        '<meta property="og:url" content="" />',
    )
    html = html.replace(
        '<title>48×48 Microwell Layout &amp; Cell Localization</title>',
        '<title>MC38 · 48×48 Microwell Localization</title>',
    )
    return html


def compare_to_published(names: np.ndarray, sample: np.ndarray, placed: dict, h5ad: str) -> None:
    adata = ad.read_h5ad(h5ad, backed="r")
    obs = adata.obs.copy()
    obs["barcode"] = adata.obs_names.astype(str)
    obs["cx"] = obs["spatial_coordinate_x"].astype(int)
    obs["cy"] = obs["spatial_coordinate_y"].astype(int)
    adata.file.close()
    df = pd.DataFrame(
        {
            "barcode": names,
            "sample_id": sample,
            "map_row": placed["map_row"].astype(int),
            "map_col": placed["map_col"].astype(int),
        }
    )
    m = df.merge(obs[["barcode", "sample", "cx", "cy"]], left_on=["barcode", "sample_id"], right_on=["barcode", "sample"])
    if m.empty:
        print("no overlap with published spatial coordinates")
        return
    exact = ((m.map_row == m.cx) & (m.map_col == m.cy)).mean()
    row = (m.map_row == m.cx).mean()
    col = (m.map_col == m.cy).mean()
    l1 = np.abs(m.map_row - m.cx) + np.abs(m.map_col - m.cy)
    shift = ((m.map_row == m.cx) & (m.map_col - 2 == m.cy)).mean()
    print(
        f"vs obsm spatial (n={len(m)} overlap): exact={exact*100:.2f}%  "
        f"row==x {row*100:.1f}%  col==y {col*100:.1f}%  "
        f"row==x & col-2==y {shift*100:.1f}%  median L1={np.median(l1):.1f}  "
        f"Pearson row-x {np.corrcoef(m.map_row, m.cx)[0,1]:.3f}  "
        f"col-y {np.corrcoef(m.map_col, m.cy)[0,1]:.3f}"
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default=str(DEFAULT_H5AD))
    ap.add_argument("--out", default=str(HERE / "sites" / "wagner_plate1"))
    ap.add_argument("--min-layout-umi", type=int, default=MIN_LAYOUT_UMI)
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    e15_bcs = np.array(gzip.open(PALAK_E15S_BCS, "rt").read().splitlines(), dtype=object)
    print(f"Palak e15s filtered GEX barcodes: {len(e15_bcs)}")
    e15_adt, adt_names, e15_found = adt_from_10x_h5(E15S_ADT_H5, e15_bcs)
    e15_bcs, e15_adt = e15_bcs[e15_found], e15_adt[e15_found]
    e15_sample = np.array(["E15S"] * len(e15_bcs), dtype=object)

    print(f"reading E14S obs from {args.h5ad}")
    adata = ad.read_h5ad(args.h5ad, backed="r")
    is_e14 = adata.obs["sample"].astype(str).to_numpy() == "E14S"
    e14_bcs = adata.obs_names.astype(str).to_numpy()[is_e14]
    adata.file.close()
    e14_adt, e14_names, e14_found = adt_from_10x_h5(E14S_ADT_H5, e14_bcs)
    if list(e14_names) != list(adt_names):
        raise SystemExit("E14S/E15S ADT feature names differ")
    e14_bcs, e14_adt = e14_bcs[e14_found], e14_adt[e14_found]
    e14_sample = np.array(["E14S"] * len(e14_bcs), dtype=object)

    names = np.concatenate([e14_bcs, e15_bcs])
    sample = np.concatenate([e14_sample, e15_sample])
    adt = np.vstack([e14_adt, e15_adt])
    print(f"combined {len(names)} cells with ADT")

    blocks = plate_blocks_from_adt(adt, adt_names, str(DEFAULT_DATA_JS), str(DEFAULT_FEATURE_REF))
    layout_umi = (
        blocks["row_p1"] + blocks["row_p2"] + blocks["col_p1"] + blocks["col_p2"]
    ).sum(axis=1)
    keep = layout_umi >= args.min_layout_umi
    blocks = {k: v[keep] for k, v in blocks.items()}
    names = names[keep]
    sample = sample[keep]
    print(f"kept {len(names)} cells with layout UMI≥{args.min_layout_umi} (both plates)")

    placed = place_from_counts(
        blocks["row_p1"],
        blocks["row_p2"],
        blocks["col_p1"],
        blocks["col_p2"],
        beta=BETA,
        use_p1=True,
        use_p2=True,
    )
    placed["layout_umi"] = (
        blocks["row_p1"] + blocks["row_p2"] + blocks["col_p1"] + blocks["col_p2"]
    ).sum(axis=1).astype(np.int64)

    ds = out_dir / "datasets" / "MC38"
    write_cells_js(ds / "cells.js", names, sample, placed, blocks)
    shutil.copyfile(DEFAULT_DATA_JS, out_dir / "data.js")

    csv_path = out_dir / "assignments.csv"
    pd.DataFrame(
        {
            "barcode": names,
            "sample_id": sample,
            "map_row": placed["map_row"].astype(int),
            "map_col": placed["map_col"].astype(int),
            "confidence": placed["confidence"],
            "row_entropy": placed["row_entropy"],
            "col_entropy": placed["col_entropy"],
            "total_entropy": placed["total_entropy"],
            "spatial_entropy": placed["spatial_entropy"],
            "layout_umi": placed["layout_umi"].astype(int),
        }
    ).to_csv(csv_path, index=False)
    print(f"wrote {csv_path}")

    html = patch_html(DEFAULT_TEMPLATE, len(names))
    (out_dir / "index.html").write_text(html)
    print(f"wrote {out_dir / 'index.html'}")
    te = placed["total_entropy"]
    print(
        f"plates=both  discrete median={np.median(te):.3f}  "
        f"spatial median={np.median(placed['spatial_entropy']):.3f}  "
        f"frac(conf>0.9)={(placed['confidence'] > 0.9).mean():.3f}  "
        f"frac(H<1)={(te < 1).mean():.3f}  "
        f"unique wells={len(set(zip(placed['map_row'], placed['map_col'])))}  "
        f"per sample {pd.Series(sample).value_counts().to_dict()}"
    )
    compare_to_published(names, sample, placed, args.h5ad)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
