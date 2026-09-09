#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PIPELINE_ROOT))

from lib.compare import compare_ocm_to_cellranger_multi, compare_to_cellranger_count
from lib.config import normalize_config
from lib.demux_ocm import demux_ocm
from lib.env import setup_env
from lib.io_utils import find_alevin_dir, find_quants_h5ad
from lib.merge import merge_gex_adt
from lib.quant import ensure_chemistry, run_simpleaf_quant, simpleaf_set_paths


def load_config(path: Path) -> dict:
    return json.loads(Path(path).read_text())


def _p(v) -> Path | None:
    return None if v in (None, "") else Path(v)


def resolve_gex_inputs(cfg: dict, outdir: Path) -> tuple[Path | None, Path | None]:
    eq = cfg.get("existing_quants") or {}
    if eq.get("gex_h5ad") and Path(eq["gex_h5ad"]).exists():
        return Path(eq["gex_h5ad"]), None
    if cfg["gex"].get("h5ad") and Path(cfg["gex"]["h5ad"]).exists():
        return Path(cfg["gex"]["h5ad"]), None
    qdir = outdir / "gex_quant"
    if (qdir / "af_quant").exists():
        return find_quants_h5ad(qdir), find_alevin_dir(qdir)
    return None, None


def resolve_adt_inputs(cfg: dict, outdir: Path) -> tuple[Path | None, Path | None]:
    if not cfg.get("adt"):
        return None, None
    eq = cfg.get("existing_quants") or {}
    if eq.get("adt_h5ad") and Path(eq["adt_h5ad"]).exists():
        return Path(eq["adt_h5ad"]), None
    if cfg["adt"].get("h5ad") and Path(cfg["adt"]["h5ad"]).exists():
        return Path(cfg["adt"]["h5ad"]), None
    qdir = outdir / "adt_quant"
    if (qdir / "af_quant").exists():
        return find_quants_h5ad(qdir), find_alevin_dir(qdir)
    return None, None


def maybe_quant(cfg: dict, outdir: Path, skip_quant: bool) -> None:
    if skip_quant:
        return
    gex = cfg["gex"]
    adt = cfg.get("adt")
    gex_ready = bool(gex.get("h5ad") and Path(gex["h5ad"]).exists())
    adt_ready = adt is None or bool(adt.get("h5ad") and Path(adt["h5ad"]).exists())
    if gex_ready and adt_ready:
        return
    threads = int(cfg.get("threads", 16))
    gex = cfg["gex"]
    chem = cfg.get("chemistries") or {}
    for name, spec in chem.items():
        ensure_chemistry(name, spec["geometry"], spec.get("expected_ori", "fw"))
    for lib in (gex, cfg.get("adt")):
        if lib and lib.get("geometry"):
            ensure_chemistry(lib["chemistry"], lib["geometry"], lib.get("expected_ori", "fw"))
    run_simpleaf_quant(
        reads1=gex["reads1"],
        reads2=gex["reads2"],
        index=Path(gex["index"]),
        chemistry=gex["chemistry"],
        output=outdir / "gex_quant",
        threads=int(gex.get("threads", threads)),
        min_reads=int(gex.get("min_reads", 10)),
        resolution=gex.get("resolution", "cr-like"),
        log_path=outdir / "logs" / "quant_gex.log",
    )
    adt = cfg.get("adt")
    if not adt:
        return
    run_simpleaf_quant(
        reads1=adt["reads1"],
        reads2=adt["reads2"],
        index=Path(adt["index"]),
        chemistry=adt["chemistry"],
        output=outdir / "adt_quant",
        threads=int(adt.get("threads", max(8, threads // 4))),
        min_reads=int(adt.get("min_reads", 10)),
        resolution=adt.get("resolution", "cr-like"),
        log_path=outdir / "logs" / "quant_adt.log",
    )


def run_pipeline(cfg: dict, *, skip_quant: bool | None = None) -> Path:
    cfg = normalize_config(cfg)
    outdir = Path(cfg["output"])
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(exist_ok=True)
    (outdir / "run_config.json").write_text(json.dumps(cfg, indent=2) + "\n")
    setup_env(Path(cfg["alevin_fry_home"]), cfg.get("tools"))
    simpleaf_set_paths()
    skip = cfg.get("skip_quant", False) if skip_quant is None else skip_quant
    maybe_quant(cfg, outdir, skip)

    gex_h5, gex_al = resolve_gex_inputs(cfg, outdir)
    adt_h5, adt_al = resolve_adt_inputs(cfg, outdir)
    if gex_h5 is None and gex_al is None:
        raise FileNotFoundError("No GEX quantification found. Provide gex.reads1/reads2/index, or gex.h5ad / existing_quants.")

    min_umi = int(cfg["min_gex_umi"])
    sample = cfg["sample"]
    feat = _p((cfg.get("adt") or {}).get("feature_ref"))
    cr = cfg.get("cellranger") or {}
    mode = cfg.get("mode", "quant")

    if mode == "ocm":
        ocm = cfg["ocm"]
        samples = ocm["samples"]
        demux_ocm(
            gex_h5ad=gex_h5,
            gex_alevin=gex_al,
            adt_h5ad=adt_h5,
            adt_alevin=adt_al,
            feature_ref=feat,
            outdir=outdir / "ocm",
            sample=sample,
            min_gex_umi=min_umi,
            samples=samples,
            overhang_start=int(ocm.get("overhang_start", 7)),
            overhang_len=int(ocm.get("overhang_len", 2)),
            include_unassigned=bool(ocm.get("include_unassigned", False)),
            cellranger_per_sample_outs=_p(cr.get("per_sample_outs")),
        )
        if cr.get("per_sample_outs"):
            compare_ocm_to_cellranger_multi(
                outdir / "ocm",
                Path(cr["per_sample_outs"]),
                outdir / "ocm" / "compare_to_cellranger.json",
                sample_ids=[s["sample_id"] for s in samples] or None,
            )
        if cr.get("filtered_mtx"):
            compare_to_cellranger_count(
                outdir / "ocm" / "union_filtered_feature_bc_matrix",
                Path(cr["filtered_mtx"]),
                outdir / "ocm" / "compare_to_cellranger_count.json",
                simpleaf_raw_barcodes=outdir / "ocm" / "union_raw_feature_bc_matrix" / "barcodes.tsv.gz",
            )
    else:
        merge_gex_adt(
            gex_h5ad=gex_h5,
            gex_alevin=gex_al,
            adt_h5ad=adt_h5,
            adt_alevin=adt_al,
            feature_ref=feat,
            outdir=outdir / "quant",
            sample=sample,
            min_gex_umi=min_umi,
            chemistry=cfg["gex"].get("chemistry", ""),
        )
        if cr.get("filtered_mtx"):
            compare_to_cellranger_count(
                outdir / "quant" / "filtered_feature_bc_matrix",
                Path(cr["filtered_mtx"]),
                outdir / "quant" / "compare_to_cellranger.json",
                simpleaf_raw_barcodes=outdir / "quant" / "raw_feature_bc_matrix" / "barcodes.tsv.gz",
            )
    return outdir


def main() -> None:
    ap = argparse.ArgumentParser(
        description="simpleaf 10x GEX (± ADT) quantification, with optional OCM demux.",
        usage="%(prog)s CONFIG.json [-o DIR] [--skip-quant]",
    )
    ap.add_argument("config", help="JSON config (see README)")
    ap.add_argument("-o", "--output", default=None, help="override config output directory")
    ap.add_argument("--skip-quant", action="store_true", help="reuse existing simpleaf quants; only merge/demux")
    args = ap.parse_args()
    cfg = load_config(Path(args.config))
    if args.output:
        cfg["output"] = args.output
    run_pipeline(cfg, skip_quant=True if args.skip_quant else None)


if __name__ == "__main__":
    main()
