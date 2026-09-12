#!/usr/bin/env python3
"""Build a here.now-ready static site from a 2-barcode layout + AnnData ADT.

Required inputs
---------------
  --h5ad          AnnData with obsm['ADT'] (and optionally obs['sample'])
  --layout        2D chip layout CSV (default: layouts/layout_2d.csv)
  --feature-ref   ADT feature reference with name,sequence columns
  --out           Output directory (index.html, data.js, datasets/<name>/cells.js,
                  assignments.csv)

Optional
--------
  --dataset       Dataset id used in the UI and cells.js path (default: cells)
  --title         Page <title> / heading
  --note          Short dataset note shown in the legend
  --min-layout-umi
  --publish       After writing files, publish with here.now (publish.sh)

Example
-------
  python build_2barcode_site.py \\
      --h5ad /path/to/gex_adt.h5ad \\
      --layout layouts/layout_2d.csv \\
      --feature-ref ../refs/new_feature_ref_quant.csv \\
      --dataset MyExp \\
      --out /tmp/myexp_site \\
      --min-layout-umi 10
"""
from __future__ import annotations

import argparse
import base64
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from cell_placement import (  # noqa: E402
    BETA,
    DEFAULT_FEATURE_REF,
    GRID,
    MAX_ENTROPY_BITS,
    peak_norm_u8,
)
from cellplacement_2barcodes import (  # noqa: E402
    DEFAULT_LAYOUT_2D,
    place_h5ad,
    sample_labels,
    write_assignments,
)

TEMPLATE = HERE / "site_template_2barcodes" / "index.html"


def write_layout_js(path: Path, barcodes: list[dict], args) -> None:
    payload = {
        "grid_size": GRID,
        "layout": "two_barcode",
        "row_plate": args.row_plate,
        "row_axis": args.row_axis,
        "col_plate": args.col_plate,
        "col_axis": args.col_axis,
        "barcodes": barcodes,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("window.LAYOUT_DATA = " + json.dumps(payload, separators=(", ", ": ")) + ";\n")
    print(f"Wrote {path} ({len(barcodes)} barcodes)")


def b64_arr(a: np.ndarray) -> str:
    return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode("ascii")


def f32_list(a: np.ndarray) -> list[float]:
    return [float(x) for x in np.asarray(a, dtype=np.float64).ravel()]


def render_index(
    template: str,
    *,
    dataset_id: str,
    dataset_label: str,
    title: str,
    heading: str,
    subtitle: str,
    note: str,
) -> str:
    repl = {
        "__SITE_TITLE__": title,
        "__SITE_HEADING__": heading,
        "__SUBTITLE__": subtitle,
        "__DATASET_ID__": dataset_id,
        "__DATASET_LABEL__": dataset_label,
        "__DATASET_NOTE__": note,
    }
    out = template
    for k, v in repl.items():
        out = out.replace(k, v)
    return out


def find_publish_sh() -> Path | None:
    env = os.environ.get("HERENOW_PUBLISH")
    if env and Path(env).is_file():
        return Path(env)
    home = Path.home()
    candidates = [
        home / ".claude/skills/here-now/scripts/publish.sh",
        home / ".cursor/skills/here-now/scripts/publish.sh",
    ]
    for p in candidates:
        if p.is_file():
            return p
    return None


def write_cells_js(path: Path, adata, placed: dict, samples: list[str] | None) -> None:
    payload = {
        "layout": "two_barcode",
        "grid_size": GRID,
        "n_cells": int(adata.n_obs),
        "max_entropy_bits": MAX_ENTROPY_BITS,
        "beta": BETA,
        "obs_names": list(map(str, adata.obs_names)),
        "map_row": placed["map_row"].astype(int).tolist(),
        "map_col": placed["map_col"].astype(int).tolist(),
        "confidence": f32_list(placed["confidence"]),
        "row_entropy": f32_list(placed["row_entropy"]),
        "col_entropy": f32_list(placed["col_entropy"]),
        "total_entropy": f32_list(placed["total_entropy"]),
        "spatial_entropy": f32_list(placed["spatial_entropy"]) if "spatial_entropy" in placed else f32_list(placed["total_entropy"]),
        "total_counts": placed["layout_umi"].astype(int).tolist(),
        "row_post_b64": b64_arr(peak_norm_u8(placed["row_post"])),
        "col_post_b64": b64_arr(peak_norm_u8(placed["col_post"])),
        "row_ll_b64": b64_arr(np.ascontiguousarray(placed["row_ll"], dtype=np.float32)),
        "col_ll_b64": b64_arr(np.ascontiguousarray(placed["col_ll"], dtype=np.float32)),
        "row_prior": f32_list(placed["row_prior"]),
        "col_prior": f32_list(placed["col_prior"]),
    }
    if samples is not None:
        payload["sample_id"] = samples
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("window.CELL_DATA = " + json.dumps(payload, separators=(",", ":")) + ";\n")
    print(f"wrote {path} n={adata.n_obs} bytes={path.stat().st_size}")


def build_site(args) -> Path:
    if not TEMPLATE.is_file():
        raise SystemExit(f"missing HTML template {TEMPLATE}")
    out = Path(args.out)
    dataset = args.dataset.strip().replace(" ", "_") or "cells"
    adata, placed, barcodes, _rows, _cols = place_h5ad(
        args.h5ad,
        args.layout,
        args.feature_ref,
        adt_names=args.adt_names,
        min_layout_umi=args.min_layout_umi,
        row_plate=args.row_plate,
        row_axis=args.row_axis,
        col_plate=args.col_plate,
        col_axis=args.col_axis,
    )
    samples = sample_labels(adata)
    n = int(adata.n_obs)
    wells = len({(int(r), int(c)) for r, c in zip(placed["map_row"], placed["map_col"])})
    title = args.title or f"{dataset} · 48×48 2-barcode localization"
    heading = args.heading or "48×48 Microwell · 2-barcode layout"
    note = args.note or "2 oligos/well · Plate-1 rows × Plate-2 columns"
    subtitle = args.subtitle or (
        f"{dataset} · {n:,} cells · {wells} occupied wells · layout UMI≥{args.min_layout_umi} · "
        f"<a href='assignments.csv' style='color:#6bcf9e'>assignments.csv</a>"
    )

    out.mkdir(parents=True, exist_ok=True)
    write_layout_js(out / "data.js", barcodes, args)
    write_cells_js(out / "datasets" / dataset / "cells.js", adata, placed, samples)
    write_assignments(out / "assignments.csv", list(map(str, adata.obs_names)), placed, samples)
    html = render_index(
        TEMPLATE.read_text(),
        dataset_id=dataset,
        dataset_label=args.dataset_label or dataset,
        title=title,
        heading=heading,
        subtitle=subtitle,
        note=note,
    )
    (out / "index.html").write_text(html)
    print(f"wrote {out / 'index.html'}")
    return out


def publish_site(site: Path, title: str, description: str) -> None:
    script = find_publish_sh()
    if script is None:
        raise SystemExit(
            "Could not find here.now publish.sh. Set HERENOW_PUBLISH=/path/to/publish.sh "
            "or pass the site directory to publish.sh yourself."
        )
    cmd = [str(script), str(site), "--client", "cursor", "--title", title, "--description", description]
    print("running", " ".join(cmd))
    subprocess.check_call(cmd)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--layout", default=str(DEFAULT_LAYOUT_2D))
    ap.add_argument("--feature-ref", default=str(DEFAULT_FEATURE_REF))
    ap.add_argument("--adt-names", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dataset", default="cells", help="Dataset id (folder name under datasets/)")
    ap.add_argument("--dataset-label", default=None, help="Label shown in the dataset dropdown")
    ap.add_argument("--title", default=None)
    ap.add_argument("--heading", default=None)
    ap.add_argument("--subtitle", default=None)
    ap.add_argument("--note", default=None)
    ap.add_argument("--min-layout-umi", type=int, default=10)
    ap.add_argument("--row-plate", type=int, default=1)
    ap.add_argument("--row-axis", default="row", choices=("row", "column"))
    ap.add_argument("--col-plate", type=int, default=2)
    ap.add_argument("--col-axis", default="column", choices=("row", "column"))
    ap.add_argument("--publish", action="store_true", help="Publish the output folder with here.now")
    args = ap.parse_args()

    site = build_site(args)
    if args.publish:
        publish_site(
            site,
            args.title or f"{args.dataset} 48x48 2-barcode localization",
            args.note or "2-barcode microwell localization viewer",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
