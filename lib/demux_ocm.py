from __future__ import annotations

import gzip
import json
from pathlib import Path

import anndata as ad
import numpy as np

from .io_utils import attach_adt, load_barcode_set, load_counts, prepare_adt, prepare_gex, strip_gem, write_mtx

GEMX_OCM_DEFAULT = [
    {"ocm_barcode_id": "OB1", "sample_id": "OB1", "description": "OCM OB1", "overhang": "GT"},
    {"ocm_barcode_id": "OB2", "sample_id": "OB2", "description": "OCM OB2", "overhang": "CA"},
    {"ocm_barcode_id": "OB3", "sample_id": "OB3", "description": "OCM OB3", "overhang": "TC"},
    {"ocm_barcode_id": "OB4", "sample_id": "OB4", "description": "OCM OB4", "overhang": "AG"},
]


def overhang_of(bc: str, start: int = 7, length: int = 2) -> str:
    b = strip_gem(bc)
    if len(b) < start + length:
        return "OTHER"
    return b[start : start + length]


def demux_ocm(
    *,
    gex_h5ad: Path | None,
    gex_alevin: Path | None,
    adt_h5ad: Path | None,
    adt_alevin: Path | None,
    feature_ref: Path | None,
    outdir: Path,
    sample: str,
    min_gex_umi: int,
    samples: list[dict],
    overhang_start: int = 7,
    overhang_len: int = 2,
    include_unassigned: bool = False,
    cellranger_per_sample_outs: Path | None = None,
) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    sample_map = list(samples or GEMX_OCM_DEFAULT)
    oh_to_ob = {s["overhang"]: s["ocm_barcode_id"] for s in sample_map}

    gex = prepare_gex(load_counts(gex_h5ad, gex_alevin, "gex"))
    gex.obs["ocm_overhang"] = [overhang_of(b, overhang_start, overhang_len) for b in gex.obs_names.astype(str)]
    gex.obs["ocm_barcode_id"] = [oh_to_ob.get(o, "OTHER") for o in gex.obs["ocm_overhang"]]

    adt = None
    if adt_h5ad is not None or adt_alevin is not None:
        try:
            adt = prepare_adt(load_counts(adt_h5ad, adt_alevin, "adt"), feature_ref)
            adt.obs["ocm_overhang"] = [overhang_of(b, overhang_start, overhang_len) for b in adt.obs_names.astype(str)]
            adt.obs["ocm_barcode_id"] = [oh_to_ob.get(o, "OTHER") for o in adt.obs["ocm_overhang"]]
        except FileNotFoundError:
            adt = None

    if include_unassigned:
        assigned = {s["ocm_barcode_id"] for s in sample_map}
        n_other = int((gex.obs["ocm_barcode_id"].values == "OTHER").sum())
        if n_other and "OTHER" not in assigned:
            sample_map.append(
                {
                    "ocm_barcode_id": "OTHER",
                    "sample_id": "unassigned",
                    "description": "barcodes without a configured OCM overhang",
                    "overhang": "*",
                }
            )

    summary = {
        "sample": sample,
        "min_gex_umi": min_gex_umi,
        "overhang_start": overhang_start,
        "overhang_len": overhang_len,
        "samples": {},
    }

    for spec in sample_map:
        ob = spec["ocm_barcode_id"]
        sample_id = spec["sample_id"]
        desc = spec.get("description", ob)
        sample_dir = outdir / "per_sample_outs" / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        gex_s = gex[gex.obs["ocm_barcode_id"].values == ob].copy()
        adt_s = None
        if adt is not None:
            adt_s = adt[adt.obs["ocm_barcode_id"].values == ob].copy()

        raw_bc_dir = sample_dir / "sample_raw_feature_bc_matrix"
        raw_bc_dir.mkdir(parents=True, exist_ok=True)
        with gzip.open(raw_bc_dir / "barcodes.tsv.gz", "wt") as fh:
            for bc in gex_s.obs_names.astype(str):
                fh.write(f"{strip_gem(bc)}-1\n")
        gex_s.write_h5ad(sample_dir / "sample_raw_gex.h5ad", compression="gzip")
        if adt_s is not None:
            adt_s.write_h5ad(sample_dir / "sample_raw_adt.h5ad", compression="gzip")

        indep = attach_adt(gex_s[gex_s.obs["gex_counts"].values >= min_gex_umi].copy(), adt_s)
        indep.obs["sample_id"] = sample_id
        indep.obs["ocm_barcode_id"] = ob
        indep.obs["description"] = desc
        indep.uns.update(
            {
                "sample": sample,
                "sample_id": sample_id,
                "ocm_barcode_id": ob,
                "filter": f"OCM {ob}; GEX UMI>={min_gex_umi}; ADT left-join",
            }
        )
        write_mtx(sample_dir / "sample_filtered_feature_bc_matrix", indep)
        indep.write_h5ad(sample_dir / "sample_filtered_gex_adt.h5ad", compression="gzip")

        cr_n = None
        cr_inter = None
        if cellranger_per_sample_outs is not None:
            cr_bc_path = (
                Path(cellranger_per_sample_outs)
                / sample_id
                / "sample_filtered_feature_bc_matrix"
                / "barcodes.tsv.gz"
            )
            if cr_bc_path.exists():
                cr_bcs = load_barcode_set(cr_bc_path)
                cr_n = len(cr_bcs)
                present = [b for b in gex_s.obs_names.astype(str) if b in cr_bcs]
                cr_inter = len(present)
                matched = attach_adt(gex_s[present].copy(), adt_s)
                matched.obs["sample_id"] = sample_id
                matched.obs["ocm_barcode_id"] = ob
                matched.obs["description"] = desc
                matched.uns.update(
                    {
                        "sample": sample,
                        "sample_id": sample_id,
                        "ocm_barcode_id": ob,
                        "filter": f"OCM {ob}; exact Cell Ranger filtered barcodes",
                    }
                )
                write_mtx(sample_dir / "sample_filtered_feature_bc_matrix_cr_match", matched)
                matched.write_h5ad(sample_dir / "sample_filtered_gex_adt_cr_match.h5ad", compression="gzip")

        summary["samples"][sample_id] = {
            "ocm_barcode_id": ob,
            "description": desc,
            "overhang": spec.get("overhang"),
            "n_raw": int(gex_s.n_obs),
            "n_filtered_independent": int(indep.n_obs),
            "n_has_adt": int(indep.obs["has_adt"].sum()) if "has_adt" in indep.obs else 0,
            "n_cellranger_filtered": cr_n,
            "n_cr_barcodes_in_simpleaf": cr_inter,
        }
        print(
            f"{sample_id} ({ob}): raw={gex_s.n_obs} "
            f"filt_umi>={min_gex_umi}={indep.n_obs} "
            f"cr_match={cr_inter}/{cr_n}"
        )

    parts = []
    for spec in sample_map:
        p = outdir / "per_sample_outs" / spec["sample_id"] / "sample_filtered_gex_adt.h5ad"
        if not p.exists():
            continue
        part = ad.read_h5ad(p)
        if part.n_obs > 0:
            parts.append(part)
    combined_path = outdir / f"{sample}_gex_adt_ocm.h5ad"
    if parts:
        combined = ad.concat(parts, join="outer", merge="same")
        combined.obs_names_make_unique()
        combined.uns["sample"] = sample
        combined.uns["ocm"] = f"overhang_pos_{overhang_start}_{overhang_start + overhang_len}"
        combined.write_h5ad(combined_path, compression="gzip")
        write_mtx(outdir / "union_filtered_feature_bc_matrix", combined)

    raw_union = outdir / "union_raw_feature_bc_matrix"
    raw_union.mkdir(parents=True, exist_ok=True)
    with gzip.open(raw_union / "barcodes.tsv.gz", "wt") as fh:
        for bc in gex.obs_names.astype(str):
            fh.write(f"{strip_gem(bc)}-1\n")

    if cellranger_per_sample_outs is not None:
        mparts = []
        for spec in sample_map:
            p = outdir / "per_sample_outs" / spec["sample_id"] / "sample_filtered_gex_adt_cr_match.h5ad"
            if p.exists():
                mparts.append(ad.read_h5ad(p))
        if mparts:
            matched = ad.concat(mparts, join="outer", merge="same")
            matched.obs_names_make_unique()
            matched.write_h5ad(outdir / f"{sample}_gex_adt_ocm_cr_match.h5ad", compression="gzip")

    (outdir / "demux_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return combined_path
