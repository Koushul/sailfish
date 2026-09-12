#!/usr/bin/env python3
"""Two-barcode 48×48 microwell placement (one row oligo × one column oligo).

Each well is the product of independent row and column posteriors. Counts on
each axis become a soft one-hot multinomial log-likelihood (β = 3) plus an
empirical log-prior, then softmax. There is no plate combining.

Default layout (`layouts/layout_2d.csv`): Plate 1 ROW 1–48 × Plate 2 COLUMN 1–48.

Examples
--------
  python cellplacement_2barcodes.py selftest

  python cellplacement_2barcodes.py from-h5ad \\
      --h5ad /path/to/gex_adt.h5ad \\
      --layout layouts/layout_2d.csv \\
      --feature-ref ../refs/new_feature_ref_quant.csv \\
      --min-layout-umi 10 \\
      --out assignments.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from cell_placement import (  # noqa: E402
    BETA,
    DEFAULT_FEATURE_REF,
    DEFAULT_SPATIAL_SIGMA,
    GRID,
    axis_loglik,
    axis_prior,
    compact_oligo_sequence,
    entropy_bits,
    extract_core,
    normalize_log_post,
    spatial_entropy_bits,
)

CHIP_LOC_RE = re.compile(r"^(ROW|COLUMN)\s+(\d+)$", re.I)
DEFAULT_LAYOUT_2D = HERE / "layouts" / "layout_2d.csv"


def place_from_axis_counts(
    row_counts: np.ndarray,
    col_counts: np.ndarray,
    beta: float = BETA,
    spatial_sigma: float | None = DEFAULT_SPATIAL_SIGMA,
) -> dict[str, np.ndarray]:
    """Place cells from one (n, 48) count matrix per axis."""
    row_counts = np.asarray(row_counts, dtype=np.float64)
    col_counts = np.asarray(col_counts, dtype=np.float64)
    if row_counts.ndim != 2 or col_counts.ndim != 2:
        raise ValueError("row_counts and col_counts must be 2-D")
    if row_counts.shape[0] != col_counts.shape[0]:
        raise ValueError("row/col count matrices must have the same number of cells")
    if row_counts.shape[1] != GRID or col_counts.shape[1] != GRID:
        raise ValueError(f"expected (n, {GRID}) count matrices")
    row_ll = axis_loglik(row_counts, beta).astype(np.float64)
    col_ll = axis_loglik(col_counts, beta).astype(np.float64)
    row_pr = axis_prior(row_counts)
    col_pr = axis_prior(col_counts)
    rp = normalize_log_post(row_ll + row_pr)
    cp = normalize_log_post(col_ll + col_pr)
    row_h = entropy_bits(rp)
    col_h = entropy_bits(cp)
    out = {
        "map_row": rp.argmax(axis=1) + 1,
        "map_col": cp.argmax(axis=1) + 1,
        "confidence": rp.max(axis=1) * cp.max(axis=1),
        "row_entropy": row_h,
        "col_entropy": col_h,
        "total_entropy": row_h + col_h,
        "row_post": rp,
        "col_post": cp,
        "row_ll": row_ll.astype(np.float32),
        "col_ll": col_ll.astype(np.float32),
        "row_prior": row_pr,
        "col_prior": col_pr,
        "layout_umi": (row_counts + col_counts).sum(axis=1).astype(np.int64),
    }
    if spatial_sigma is not None:
        out["spatial_entropy"] = spatial_entropy_bits(rp, cp, sigma=spatial_sigma)
        out["spatial_sigma"] = np.asarray(spatial_sigma, dtype=np.float64)
    return out


def layout_barcodes_from_csv(path: str | Path) -> list[dict]:
    """Parse a 2D chip layout CSV into site `data.js` barcode records."""
    try:
        import pandas as pd
    except ImportError as e:
        raise SystemExit(f"layout CSV requires pandas: {e}") from e

    df = pd.read_csv(path)
    needed = ["Name", "Chip Location (Row/Column)"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise SystemExit(f"{path} missing columns {missing}; have {list(df.columns)}")
    df = df.dropna(subset=needed)
    barcodes = []
    for _, row in df.iterrows():
        loc = str(row["Chip Location (Row/Column)"]).strip()
        m = CHIP_LOC_RE.match(loc)
        if not m:
            raise SystemExit(f"Unparseable chip location {loc!r} for {row['Name']}")
        plate_raw = row["96 Source Plate #"] if "96 Source Plate #" in df.columns else 1
        plate_s = str(plate_raw).replace("Plate", "").strip()
        well = (
            str(row["Well Position in 384 Source Plate"])
            if "Well Position in 384 Source Plate" in df.columns
            else ""
        )
        barcodes.append(
            {
                "name": str(row["Name"]),
                "plate": int(plate_s),
                "axis": "row" if m.group(1).upper() == "ROW" else "column",
                "index": int(m.group(2)),
                "well": well,
                "sequence": compact_oligo_sequence(row["Sequence"]) if "Sequence" in df.columns else "",
            }
        )
    return barcodes


def two_barcode_families(
    barcodes: list[dict],
    row_plate: int = 1,
    row_axis: str = "row",
    col_plate: int = 2,
    col_axis: str = "column",
) -> tuple[list[dict], list[dict]]:
    """Select the 48 row oligos and 48 column oligos that define each well."""
    row_axis = "row" if row_axis in ("row", "ROW") else "column"
    col_axis = "row" if col_axis in ("row", "ROW") else "column"
    rows: dict[int, dict] = {}
    cols: dict[int, dict] = {}
    for b in barcodes:
        item = {**b, "core": extract_core(b["sequence"]) if b.get("sequence") else ""}
        if b["plate"] == row_plate and b["axis"] == row_axis:
            rows[int(b["index"])] = item
        elif b["plate"] == col_plate and b["axis"] == col_axis:
            cols[int(b["index"])] = item
    missing_r = [i for i in range(1, GRID + 1) if i not in rows]
    missing_c = [i for i in range(1, GRID + 1) if i not in cols]
    if missing_r or missing_c:
        raise SystemExit(
            f"layout missing {row_axis} (plate {row_plate}) indices {missing_r} "
            f"or {col_axis} (plate {col_plate}) indices {missing_c}"
        )
    return [rows[i] for i in range(1, GRID + 1)], [cols[i] for i in range(1, GRID + 1)]


def load_two_barcode_layout(
    path: str | Path,
    row_plate: int = 1,
    row_axis: str = "row",
    col_plate: int = 2,
    col_axis: str = "column",
) -> tuple[list[dict], list[dict], list[dict]]:
    barcodes = layout_barcodes_from_csv(path)
    rows, cols = two_barcode_families(barcodes, row_plate, row_axis, col_plate, col_axis)
    return barcodes, rows, cols


def adt_axis_counts(
    adt: np.ndarray,
    adt_names: list[str],
    barcodes: list[dict],
    feature_ref: str | Path,
) -> np.ndarray:
    """Pull (n, 48) counts for a barcode family by matching oligo 15-mers to ADT features."""
    try:
        import pandas as pd
    except ImportError as e:
        raise SystemExit(f"from-h5ad requires pandas: {e}") from e

    ref = pd.read_csv(feature_ref)
    if "name" not in ref.columns or "sequence" not in ref.columns:
        raise SystemExit(f"feature ref {feature_ref} needs 'name' and 'sequence' columns")
    name_to_seq = dict(zip(ref["name"].astype(str), ref["sequence"].astype(str)))
    seq_to_col: dict[str, int] = {}
    name_to_col: dict[str, int] = {}
    for j, n in enumerate(adt_names):
        name_to_col[str(n)] = j
        seq = name_to_seq.get(str(n))
        if seq:
            seq_to_col[seq] = j
            seq_to_col[seq.upper()] = j
    cols = []
    for bc in barcodes:
        core = bc.get("core") or (extract_core(bc["sequence"]) if bc.get("sequence") else "")
        j = seq_to_col.get(core) or seq_to_col.get(core.upper())
        if j is None:
            j = name_to_col.get(bc["name"])
        if j is None:
            raise SystemExit(f"No ADT feature for {bc['name']} core={core!r}")
        cols.append(j)
    return np.asarray(adt, dtype=np.float64)[:, cols]


def _adt_matrix_and_names(adata, feature_ref: str | Path, adt_names_path: str | None):
    from scipy import sparse

    if "ADT" not in adata.obsm:
        raise SystemExit("h5ad missing obsm['ADT']")
    adt = adata.obsm["ADT"]
    adt = adt.toarray() if sparse.issparse(adt) else np.asarray(adt, dtype=np.float64)
    if adt_names_path:
        names = Path(adt_names_path).read_text().splitlines()
        if len(names) != adt.shape[1]:
            raise SystemExit(f"--adt-names has {len(names)} lines, ADT has {adt.shape[1]} columns")
        return adt, names
    av = adata.uns.get("ADT_var")
    if av is not None:
        import pandas as pd

        df = pd.DataFrame(av) if isinstance(av, dict) else av
        col = "feature_name" if "feature_name" in df.columns else df.columns[0]
        names = df[col].astype(str).tolist()
        if len(names) != adt.shape[1]:
            raise SystemExit(f"uns['ADT_var'] length {len(names)} != ADT columns {adt.shape[1]}")
        return adt, names
    import pandas as pd

    ref = pd.read_csv(feature_ref)
    if len(ref) != adt.shape[1]:
        raise SystemExit(
            f"h5ad has no uns['ADT_var'] and feature ref has {len(ref)} rows vs ADT {adt.shape[1]} columns; "
            "pass --adt-names"
        )
    return adt, ref["name"].astype(str).tolist()


def sample_labels(adata) -> list[str] | None:
    for key in ("sample", "sample_id", "batch"):
        if key in adata.obs.columns:
            return adata.obs[key].astype(str).tolist()
    return None


def place_h5ad(
    h5ad_path: str | Path,
    layout_csv: str | Path,
    feature_ref: str | Path,
    *,
    adt_names: str | None = None,
    min_layout_umi: int = 0,
    beta: float = BETA,
    spatial_sigma: float | None = DEFAULT_SPATIAL_SIGMA,
    row_plate: int = 1,
    row_axis: str = "row",
    col_plate: int = 2,
    col_axis: str = "column",
) -> tuple[object, dict[str, np.ndarray], list[dict], list[dict], list[dict]]:
    try:
        import anndata as ad
    except ImportError as e:
        raise SystemExit(f"from-h5ad requires anndata: {e}") from e

    barcodes, row_bc, col_bc = load_two_barcode_layout(
        layout_csv, row_plate, row_axis, col_plate, col_axis
    )
    adata = ad.read_h5ad(h5ad_path)
    adt, names = _adt_matrix_and_names(adata, feature_ref, adt_names)
    row_counts = adt_axis_counts(adt, names, row_bc, feature_ref)
    col_counts = adt_axis_counts(adt, names, col_bc, feature_ref)
    placed = place_from_axis_counts(row_counts, col_counts, beta=beta, spatial_sigma=spatial_sigma)
    keep = placed["layout_umi"] >= min_layout_umi
    if not np.all(keep):
        n0 = int(keep.shape[0])
        adata = adata[keep].copy()
        placed = {
            k: (v[keep] if isinstance(v, np.ndarray) and v.ndim >= 1 and v.shape[0] == n0 else v)
            for k, v in placed.items()
        }
    return adata, placed, barcodes, row_bc, col_bc


def write_assignments(
    path: Path,
    names: list[str],
    placed: dict[str, np.ndarray],
    samples: list[str] | None = None,
) -> None:
    fields = [
        "barcode",
        "map_row",
        "map_col",
        "confidence",
        "row_entropy",
        "col_entropy",
        "total_entropy",
    ]
    if samples is not None:
        fields.insert(1, "sample_id")
    if "spatial_entropy" in placed:
        fields.append("spatial_entropy")
    if "layout_umi" in placed:
        fields.append("layout_umi")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for i, name in enumerate(names):
            row = {
                "barcode": name,
                "map_row": int(placed["map_row"][i]),
                "map_col": int(placed["map_col"][i]),
                "confidence": float(placed["confidence"][i]),
                "row_entropy": float(placed["row_entropy"][i]),
                "col_entropy": float(placed["col_entropy"][i]),
                "total_entropy": float(placed["total_entropy"][i]),
            }
            if samples is not None:
                row["sample_id"] = samples[i]
            if "spatial_entropy" in placed:
                row["spatial_entropy"] = float(placed["spatial_entropy"][i])
            if "layout_umi" in placed:
                row["layout_umi"] = int(placed["layout_umi"][i])
            w.writerow(row)
    print(f"Wrote {path} (n={len(names)})")


def _selftest() -> None:
    rng = np.random.default_rng(0)
    n = 32
    row = np.zeros((n, GRID))
    col = np.zeros((n, GRID))
    true_r = rng.integers(0, GRID, size=n)
    true_c = rng.integers(0, GRID, size=n)
    for i in range(n):
        row[i, true_r[i]] = 40
        col[i, true_c[i]] = 40
        row[i] += rng.integers(0, 2, GRID)
        col[i] += rng.integers(0, 2, GRID)
    out = place_from_axis_counts(row, col, beta=3.0)
    assert np.all(out["map_row"] == true_r + 1)
    assert np.all(out["map_col"] == true_c + 1)
    assert np.all(out["total_entropy"] < 0.5)
    assert np.all(out["confidence"] > 0.9)
    uni = np.full(GRID, 1.0 / GRID)
    assert abs(float(entropy_bits(uni[None, :])[0]) - math.log2(GRID)) < 1e-12
    csv_path = DEFAULT_LAYOUT_2D
    if csv_path.exists():
        barcodes, rows, cols = load_two_barcode_layout(csv_path)
        assert len(rows) == GRID and len(cols) == GRID
        assert len(barcodes) >= 2 * GRID
        print(f"selftest: ok  layout={csv_path.name} n_oligos={len(barcodes)}")
    else:
        print("selftest: ok  (layout_2d.csv not present, skipped CSV parse)")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    sub.add_parser("selftest", help="Synthetic 2-barcode placement checks")

    p_l = sub.add_parser("dump-layout", help="Print the 48×48 row/column oligo families from a layout CSV")
    p_l.add_argument("--layout", default=str(DEFAULT_LAYOUT_2D))
    p_l.add_argument("--row-plate", type=int, default=1)
    p_l.add_argument("--row-axis", default="row", choices=("row", "column"))
    p_l.add_argument("--col-plate", type=int, default=2)
    p_l.add_argument("--col-axis", default="column", choices=("row", "column"))

    p_a = sub.add_parser("from-counts", help="Place cells from (n, 48) .npy count matrices")
    p_a.add_argument("--row-npy", required=True)
    p_a.add_argument("--col-npy", required=True)
    p_a.add_argument("--names", default=None, help="Text file of cell barcodes (one per line)")
    p_a.add_argument("--beta", type=float, default=BETA)
    p_a.add_argument("--spatial-sigma", type=float, default=DEFAULT_SPATIAL_SIGMA)
    p_a.add_argument("--out", required=True)

    p_h = sub.add_parser("from-h5ad", help="Place cells from AnnData obsm['ADT'] + layout_2d.csv")
    p_h.add_argument("--h5ad", required=True)
    p_h.add_argument("--layout", default=str(DEFAULT_LAYOUT_2D))
    p_h.add_argument("--feature-ref", default=str(DEFAULT_FEATURE_REF))
    p_h.add_argument("--adt-names", default=None)
    p_h.add_argument("--beta", type=float, default=BETA)
    p_h.add_argument("--spatial-sigma", type=float, default=DEFAULT_SPATIAL_SIGMA)
    p_h.add_argument("--min-layout-umi", type=int, default=0)
    p_h.add_argument("--row-plate", type=int, default=1)
    p_h.add_argument("--row-axis", default="row", choices=("row", "column"))
    p_h.add_argument("--col-plate", type=int, default=2)
    p_h.add_argument("--col-axis", default="column", choices=("row", "column"))
    p_h.add_argument("--out", required=True)

    args = ap.parse_args(argv)

    if args.cmd == "selftest":
        _selftest()
        return 0

    if args.cmd == "dump-layout":
        barcodes, rows, cols = load_two_barcode_layout(
            args.layout, args.row_plate, args.row_axis, args.col_plate, args.col_axis
        )
        print(
            json.dumps(
                {
                    "n_oligos": len(barcodes),
                    "row": [{"index": b["index"], "name": b["name"], "plate": b["plate"]} for b in rows],
                    "col": [{"index": b["index"], "name": b["name"], "plate": b["plate"]} for b in cols],
                },
                indent=2,
            )
        )
        return 0

    if args.cmd == "from-counts":
        row = np.load(args.row_npy)
        col = np.load(args.col_npy)
        placed = place_from_axis_counts(row, col, beta=args.beta, spatial_sigma=args.spatial_sigma)
        if args.names:
            names = Path(args.names).read_text().splitlines()
        else:
            names = [f"cell{i}" for i in range(row.shape[0])]
        write_assignments(Path(args.out), names, placed)
        return 0

    adata, placed, _barcodes, _rows, _cols = place_h5ad(
        args.h5ad,
        args.layout,
        args.feature_ref,
        adt_names=args.adt_names,
        min_layout_umi=args.min_layout_umi,
        beta=args.beta,
        spatial_sigma=args.spatial_sigma,
        row_plate=args.row_plate,
        row_axis=args.row_axis,
        col_plate=args.col_plate,
        col_axis=args.col_axis,
    )
    write_assignments(Path(args.out), list(map(str, adata.obs_names)), placed, sample_labels(adata))
    te = placed["total_entropy"]
    wells = len({(int(r), int(c)) for r, c in zip(placed["map_row"], placed["map_col"])})
    print(
        f"n={len(te)} wells={wells} median_H={float(np.median(te)):.3f} "
        f"frac(conf>0.9)={float(np.mean(placed['confidence'] > 0.9)):.3f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
