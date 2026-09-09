#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PIPELINE_ROOT))

from lib.compare import compare_ocm_to_cellranger_multi, compare_to_cellranger_count
from lib.demux_ocm import GEMX_OCM_DEFAULT, demux_ocm
from lib.env import setup_env
from lib.io_utils import find_alevin_dir, find_quants_h5ad
from lib.merge import merge_gex_adt
from lib.quant import ensure_chemistry, run_simpleaf_quant, simpleaf_set_paths


def load_config(path: Path) -> dict:
    return json.loads(Path(path).read_text())


def _p(v) -> Path | None:
    return None if v is None else Path(v)


def resolve_gex_inputs(cfg: dict, outdir: Path) -> tuple[Path | None, Path | None]:
    eq = cfg.get("existing_quants") or {}
    if eq.get("gex_h5ad") and Path(eq["gex_h5ad"]).exists():
        return Path(eq["gex_h5ad"]), None
    if eq.get("gex_alevin") and Path(eq["gex_alevin"]).exists():
        return None, Path(eq["gex_alevin"])
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
    if eq.get("adt_alevin") and Path(eq["adt_alevin"]).exists():
        return None, Path(eq["adt_alevin"])
    qdir = outdir / "adt_quant"
    if (qdir / "af_quant").exists():
        return find_quants_h5ad(qdir), find_alevin_dir(qdir)
    return None, None


def maybe_quant(cfg: dict, outdir: Path, skip_quant: bool) -> None:
    if skip_quant:
        return
    threads = int(cfg.get("threads", 16))
    gex = cfg["gex"]
    t_gex = int(gex.get("threads", threads))
    chem = cfg.get("chemistries") or {}
    for name, spec in chem.items():
        ensure_chemistry(name, spec["geometry"], spec.get("expected_ori", "fw"))
    run_simpleaf_quant(
        reads1=gex["reads1"],
        reads2=gex["reads2"],
        index=Path(gex["index"]),
        chemistry=gex["chemistry"],
        output=outdir / "gex_quant",
        threads=t_gex,
        min_reads=int(gex.get("min_reads", 10)),
        resolution=gex.get("resolution", "cr-like"),
        log_path=outdir / "logs" / "quant_gex.log",
    )
    adt = cfg.get("adt")
    if not adt:
        return
    t_adt = int(adt.get("threads", max(8, threads // 4)))
    run_simpleaf_quant(
        reads1=adt["reads1"],
        reads2=adt["reads2"],
        index=Path(adt["index"]),
        chemistry=adt["chemistry"],
        output=outdir / "adt_quant",
        threads=t_adt,
        min_reads=int(adt.get("min_reads", 10)),
        resolution=adt.get("resolution", "cr-like"),
        log_path=outdir / "logs" / "quant_adt.log",
    )


def run_pipeline(cfg: dict, *, skip_quant: bool | None = None) -> Path:
    outdir = Path(cfg["output"])
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(exist_ok=True)
    (outdir / "run_config.json").write_text(json.dumps(cfg, indent=2) + "\n")
    setup_env(Path(cfg["alevin_fry_home"]))
    simpleaf_set_paths()
    skip = cfg.get("skip_quant", False) if skip_quant is None else skip_quant
    maybe_quant(cfg, outdir, skip)

    gex_h5, gex_al = resolve_gex_inputs(cfg, outdir)
    adt_h5, adt_al = resolve_adt_inputs(cfg, outdir)
    if gex_h5 is None and gex_al is None:
        raise FileNotFoundError("No GEX quantification found; run without --skip-quant or set existing_quants")

    filt = cfg.get("filter") or {}
    min_umi = int(filt.get("min_gex_umi", 500))
    sample = cfg.get("sample", "sample")
    feat = None
    if cfg.get("adt"):
        feat = _p(cfg["adt"].get("feature_ref"))
    mode = cfg.get("mode", "quant")

    if mode == "ocm":
        ocm = cfg.get("ocm") or {}
        cr = cfg.get("cellranger") or {}
        demux_ocm(
            gex_h5ad=gex_h5,
            gex_alevin=gex_al,
            adt_h5ad=adt_h5,
            adt_alevin=adt_al,
            feature_ref=feat,
            outdir=outdir / "ocm",
            sample=sample,
            min_gex_umi=min_umi,
            samples=ocm.get("samples") or GEMX_OCM_DEFAULT,
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
                sample_ids=[s["sample_id"] for s in (ocm.get("samples") or [])] or None,
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
        cr = cfg.get("cellranger") or {}
        if cr.get("filtered_mtx"):
            compare_to_cellranger_count(
                outdir / "quant" / "filtered_feature_bc_matrix",
                Path(cr["filtered_mtx"]),
                outdir / "quant" / "compare_to_cellranger.json",
                simpleaf_raw_barcodes=outdir / "quant" / "raw_feature_bc_matrix" / "barcodes.tsv.gz",
            )
    return outdir


def run_cellranger_count(cfg: dict) -> Path:
    cr = cfg["cellranger"]
    cellranger = cr.get("bin", "/software/rhel9/manual/install/cellranger/cellranger-10.0.0/cellranger")
    workdir = Path(cr.get("workdir", cfg["output"]))
    workdir.mkdir(parents=True, exist_ok=True)
    run_id = cr.get("id", "cr_count")
    cmd = [
        cellranger,
        "count",
        f"--id={run_id}",
        f"--transcriptome={cr['transcriptome']}",
        f"--fastqs={cr['fastqs']}",
        f"--localcores={cr.get('localcores', cfg.get('threads', 16))}",
        f"--localmem={cr.get('localmem', 64)}",
        "--create-bam=false",
    ]
    if cr.get("sample"):
        cmd.append(f"--sample={cr['sample']}")
    env = os.environ.copy()
    env.pop("PYTHONNOUSERSITE", None)
    subprocess.run(cmd, cwd=workdir, check=True, env=env)
    return workdir / run_id / "outs" / "filtered_feature_bc_matrix"


def main() -> None:
    ap = argparse.ArgumentParser(description="Dataset-agnostic simpleaf 10x GEX±ADT pipeline (quant or OCM).")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_run = sub.add_parser("run", help="Run quant or OCM workflow from a JSON config")
    p_run.add_argument("--config", required=True)
    p_run.add_argument("--skip-quant", action="store_true")
    p_run.add_argument("--output", default=None)

    p_cmp = sub.add_parser("compare-count", help="Compare simpleaf MTX to cellranger count MTX")
    p_cmp.add_argument("--simpleaf-mtx", required=True)
    p_cmp.add_argument("--cellranger-mtx", required=True)
    p_cmp.add_argument("--out-json", required=True)
    p_cmp.add_argument("--simpleaf-raw-barcodes", default=None)

    p_cr = sub.add_parser("cellranger-count", help="Run cellranger count using paths in the config")
    p_cr.add_argument("--config", required=True)

    args = ap.parse_args()
    if args.cmd == "run":
        cfg = load_config(Path(args.config))
        if args.output:
            cfg["output"] = args.output
        run_pipeline(cfg, skip_quant=True if args.skip_quant else None)
    elif args.cmd == "compare-count":
        compare_to_cellranger_count(
            Path(args.simpleaf_mtx),
            Path(args.cellranger_mtx),
            Path(args.out_json),
            Path(args.simpleaf_raw_barcodes) if args.simpleaf_raw_barcodes else None,
        )
    elif args.cmd == "cellranger-count":
        run_cellranger_count(load_config(Path(args.config)))


if __name__ == "__main__":
    main()
