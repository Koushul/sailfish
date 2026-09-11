#!/usr/bin/env python3
"""Place WagnerCollab MC38 cells with Plate 1 only and write a lucid-crystal-style site."""
from __future__ import annotations

import argparse
import base64
import gzip
import json
import shutil
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse

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
DEFAULT_CR = {
    "E14S": Path("/ix1/ylee/shared/MC38_Hypoxia_001/E14S/filtered_feature_bc_matrix"),
    "E15S": Path("/ix1/ylee/shared/MC38_Hypoxia_001/E15S/filtered_feature_bc_matrix"),
}
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


def load_10x_adt(mtx_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    barcodes = np.array(
        gzip.open(mtx_dir / "barcodes.tsv.gz", "rt").read().splitlines(), dtype=object
    )
    feats = []
    with gzip.open(mtx_dir / "features.tsv.gz", "rt") as fh:
        for line in fh:
            fid, name, ftype = line.rstrip("\n").split("\t")
            feats.append((fid, name, ftype))
    mat = mmread(mtx_dir / "matrix.mtx.gz").tocsr()
    adt_idx = [i for i, row in enumerate(feats) if row[2] == "Antibody Capture"]
    if not adt_idx:
        raise SystemExit(f"No Antibody Capture features in {mtx_dir}")
    names = np.array([feats[i][1] for i in adt_idx], dtype=object)
    adt = mat[adt_idx, :].T.tocsr()
    return barcodes, names, adt


def adt_for_velocity_cells(
    adata: ad.AnnData, cr_dirs: dict[str, Path]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sample = adata.obs["sample"].astype(str).to_numpy()
    obs_names = adata.obs_names.astype(str).to_numpy()
    n = adata.n_obs
    adt_names = None
    blocks = {}
    for lab, path in cr_dirs.items():
        bc, names, mat = load_10x_adt(path)
        if adt_names is None:
            adt_names = names
        elif list(names) != list(adt_names):
            raise SystemExit(f"ADT feature names differ for {lab}")
        lookup = {b: i for i, b in enumerate(bc)}
        blocks[lab] = (lookup, mat)

    n_adt = len(adt_names)
    out = np.zeros((n, n_adt), dtype=np.float64)
    found = np.zeros(n, dtype=bool)
    for i, (name, lab) in enumerate(zip(obs_names, sample)):
        lookup, mat = blocks[lab]
        j = lookup.get(name)
        if j is None:
            j = lookup.get(name.split("-")[0] + "-1")
        if j is None:
            continue
        row = mat[j]
        out[i] = row.toarray().ravel() if sparse.issparse(row) else np.asarray(row).ravel()
        found[i] = True
    print(f"ADT matched {int(found.sum())} / {n} velocity cells")
    return out, adt_names, found


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
          <option value="MC38" selected>MC38 Plate 1</option>
        </select>""",
    )
    html = html.replace(
        """    const DATASETS = {
      E28S: { path: 'datasets/E28S/cells.js', note: 'edge/core hypoxia OCM · original published assignment' },
      E14S: { path: 'datasets/E14S/cells.js', note: 'entropy assignment · layout UMI≥10, GEX≥500' },
      E15S: { path: 'datasets/E15S/cells.js', note: 'entropy assignment · layout UMI≥10, GEX≥500' },
    };""",
        """    const DATASETS = {
      MC38: { path: 'datasets/MC38/cells.js', note: 'WagnerCollab mc38_velocity · Plate 1 only · layout UMI≥10' },
    };""",
    )
    html = html.replace("let currentDataset = 'E28S';", "let currentDataset = 'MC38';")
    html = html.replace(
        "Hover a microwell to inspect barcodes; switch to Cells to position all cells and assess entropy",
        "WagnerCollab MC38 (E14S+E15S) · Plate 1 localization · "
        f"{n_cells:,} cells (layout UMI≥{MIN_LAYOUT_UMI}) · "
        "<a href='assignments.csv' style='color:#6bcf9e'>assignments.csv</a>",
    )
    html = html.replace(
        '<label class="checkline"><input type="checkbox" id="plate-1" checked /><span class="dot p1"></span> Plate 1</label>\n'
        '        <label class="checkline"><input type="checkbox" id="plate-2" checked /><span class="dot p2"></span> Plate 2</label>',
        '<label class="checkline"><input type="checkbox" id="plate-1" checked disabled /><span class="dot p1"></span> Plate 1</label>\n'
        '        <label class="checkline"><input type="checkbox" id="plate-2" disabled /><span class="dot p2"></span> Plate 2 (not used)</label>',
    )
    html = html.replace("plate2: true,", "plate2: false,")
    html = html.replace(
        "Plate 1 and Plate 2 each contribute independent row/column sets (192 oligos total).",
        "Localization uses Plate 1 row/column spatial-hash oligos only (96 oligos).",
    )
    html = html.replace(
        '<meta property="og:url" content="https://lucid-crystal-kmqy.here.now/" />',
        '<meta property="og:url" content="" />',
    )
    html = html.replace(
        '<title>48×48 Microwell Layout &amp; Cell Localization</title>',
        '<title>MC38 Plate 1 · 48×48 Microwell Localization</title>',
    )
    return html


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default=str(DEFAULT_H5AD))
    ap.add_argument("--out", default=str(HERE / "sites" / "wagner_plate1"))
    ap.add_argument("--min-layout-umi", type=int, default=MIN_LAYOUT_UMI)
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"reading obs from {args.h5ad}")
    adata = ad.read_h5ad(args.h5ad, backed="r")
    adt, adt_names, found = adt_for_velocity_cells(adata, DEFAULT_CR)
    names = adata.obs_names.astype(str).to_numpy()
    sample = adata.obs["sample"].astype(str).to_numpy()
    adata.file.close()

    adt = adt[found]
    names = names[found]
    sample = sample[found]
    blocks = plate_blocks_from_adt(adt, adt_names, str(DEFAULT_DATA_JS), str(DEFAULT_FEATURE_REF))
    layout_umi = (blocks["row_p1"] + blocks["col_p1"]).sum(axis=1)
    keep = layout_umi >= args.min_layout_umi
    blocks = {k: v[keep] for k, v in blocks.items()}
    names = names[keep]
    sample = sample[keep]
    print(f"kept {len(names)} cells with Plate-1 layout UMI≥{args.min_layout_umi}")

    placed = place_from_counts(
        blocks["row_p1"],
        blocks["row_p2"],
        blocks["col_p1"],
        blocks["col_p2"],
        beta=BETA,
        use_p1=True,
        use_p2=False,
    )
    placed["layout_umi"] = (blocks["row_p1"] + blocks["col_p1"]).sum(axis=1).astype(np.int64)

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
        f"plates=1  discrete median={np.median(te):.3f}  "
        f"spatial median={np.median(placed['spatial_entropy']):.3f}  "
        f"frac(conf>0.9)={(placed['confidence'] > 0.9).mean():.3f}  "
        f"frac(H<1)={(te < 1).mean():.3f}  "
        f"per sample {pd.Series(sample).value_counts().to_dict()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
