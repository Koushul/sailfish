from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from .io_utils import find_alevin_dir, find_quants_h5ad, load_counts, strip_gem
from .quant import run_simpleaf_quant
from .vdj_ref import ensure_vdj_index

TCR_LOCI = ("TRA", "TRB", "TRG", "TRD")
BCR_LOCI = ("IGH", "IGK", "IGL")


def _as_list(v) -> list[str]:
    if v is None:
        return []
    if isinstance(v, str):
        return [v]
    return [str(p) for p in v]


def resolve_vdj_fastqs(lib: dict | None) -> dict | None:
    """Return lib with reads1/reads2 filled from an Illumina FASTQ directory if needed."""
    if not lib:
        return None
    lib = dict(lib)
    r1 = _as_list(lib.get("reads1"))
    r2 = _as_list(lib.get("reads2"))
    if r1 and r2:
        lib["reads1"], lib["reads2"] = r1, r2
        return lib
    fastq_dir = lib.get("fastqs")
    if not fastq_dir:
        return lib
    d = Path(fastq_dir)
    if not d.is_dir():
        return lib
    r1 = sorted(str(p) for p in d.glob("*_R1_001.fastq.gz"))
    r2 = sorted(str(p) for p in d.glob("*_R2_001.fastq.gz"))
    if not r1:
        r1 = sorted(str(p) for p in d.glob("*_R1*.fastq.gz")) + sorted(str(p) for p in d.glob("*_R1*.fq.gz"))
        r2 = sorted(str(p) for p in d.glob("*_R2*.fastq.gz")) + sorted(str(p) for p in d.glob("*_R2*.fq.gz"))
    lib["reads1"], lib["reads2"] = r1, r2
    return lib


def _reads_exist(lib: dict | None) -> bool:
    lib = resolve_vdj_fastqs(lib)
    if not lib:
        return False
    r1, r2 = lib.get("reads1") or [], lib.get("reads2") or []
    if not r1 or not r2 or len(r1) != len(r2):
        return False
    return all(Path(p).exists() for p in list(r1) + list(r2))


def _segment_lookup(segments_tsv: Path) -> dict[str, dict]:
    df = pd.read_csv(segments_tsv, sep="\t")
    by: dict[str, dict] = {}
    for rec in df.to_dict("records"):
        by[str(rec["feature_id"])] = rec
        by[str(rec["gene_name"])] = rec
    return by


def _region_kind(region_type: str) -> str:
    if "V-REGION" in region_type:
        return "V"
    if region_type == "J-REGION":
        return "J"
    if region_type == "C-REGION":
        return "C"
    return "other"


def annotate_from_h5ad(adata: ad.AnnData, segments_tsv: Path, receptor: str) -> pd.DataFrame:
    lookup = _segment_lookup(segments_tsv)
    genes = adata.var_names.astype(str)
    kinds = []
    loci = []
    for g in genes:
        rec = lookup.get(g)
        if rec is not None:
            kinds.append(_region_kind(str(rec["region_type"])))
            loci.append(str(rec["locus"]))
        else:
            kinds.append("other")
            loci.append("")
    adata.var["vdj_kind"] = kinds
    adata.var["vdj_locus"] = loci
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X)
    rows = []
    bcs = [strip_gem(b) for b in adata.obs_names.astype(str)]
    loci_keep = TCR_LOCI if receptor == "TCR" else BCR_LOCI
    for i, bc in enumerate(bcs):
        rec = {"barcode": bc}
        total = 0.0
        for loc in loci_keep:
            mask = np.array([l == loc for l in loci])
            if not mask.any():
                rec[f"{loc}_v"] = ""
                rec[f"{loc}_j"] = ""
                rec[f"{loc}_c"] = ""
                rec[f"{loc}_umi"] = 0
                rec[f"{loc}_v_umi"] = 0
                rec[f"{loc}_j_umi"] = 0
                continue
            umi_loc = float(X[i, mask].sum())
            total += umi_loc
            rec[f"{loc}_umi"] = umi_loc
            for kind, key in (("V", "v"), ("J", "j"), ("C", "c")):
                km = mask & (np.array(kinds) == kind)
                if not km.any() or X[i, km].sum() <= 0:
                    rec[f"{loc}_{key}"] = ""
                    rec[f"{loc}_{key}_umi"] = 0
                    continue
                j = int(np.argmax(np.where(km, X[i], -1)))
                rec[f"{loc}_{key}"] = genes[j]
                rec[f"{loc}_{key}_umi"] = float(X[i, j])
        rec["total_umi"] = total
        if receptor == "TCR":
            rec["paired"] = bool(rec.get("TRA_v") and rec.get("TRA_j") and rec.get("TRB_v") and rec.get("TRB_j"))
        else:
            rec["paired"] = bool(
                rec.get("IGH_v")
                and rec.get("IGH_j")
                and (
                    (rec.get("IGK_v") and rec.get("IGK_j"))
                    or (rec.get("IGL_v") and rec.get("IGL_j"))
                )
            )
        rows.append(rec)
    return pd.DataFrame(rows)


def attach_vdj_obs(adata: ad.AnnData, ann: pd.DataFrame, prefix: str) -> ad.AnnData:
    ann = ann.copy()
    ann.index = ann["barcode"].map(strip_gem)
    map_bc = pd.Index([strip_gem(b) for b in adata.obs_names.astype(str)])
    for col in ann.columns:
        if col == "barcode":
            continue
        adata.obs[f"{prefix}_{col}"] = ann.reindex(map_bc)[col].to_numpy()
    return adata


def _reload_library(output: Path, receptor: str, segments_tsv: Path) -> dict | None:
    if not (output / "af_quant").exists():
        return None
    summary_path = output / "vdj_summary.json"
    csv_path = output / f"{receptor.lower()}_annotations.csv"
    if summary_path.exists() and csv_path.exists():
        return json.loads(summary_path.read_text())
    h5 = find_quants_h5ad(output)
    adata = load_counts(h5, find_alevin_dir(output), f"vdj_{receptor}")
    ann = annotate_from_h5ad(adata, segments_tsv, receptor)
    csv_path = output / f"{receptor.lower()}_annotations.csv"
    ann.to_csv(csv_path, index=False)
    summary = {
        "receptor": receptor,
        "n_barcodes": int(adata.n_obs),
        "n_features": int(adata.n_vars),
        "n_paired": int(ann["paired"].sum()) if "paired" in ann.columns else 0,
        "median_umi": float(np.median(ann["total_umi"])) if len(ann) else 0.0,
        "annotations": str(csv_path),
        "h5ad": str(h5),
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def run_vdj_library(
    *,
    lib: dict,
    receptor: str,
    index: Path,
    t2g: Path,
    segments_tsv: Path,
    chemistry: str,
    output: Path,
    threads: int,
    min_reads: int,
    resolution: str,
) -> dict:
    lib = resolve_vdj_fastqs(lib)
    existing = _reload_library(output, receptor, segments_tsv)
    if existing:
        print(f"vdj: reusing {receptor} quant at {output}")
        return existing
    output.mkdir(parents=True, exist_ok=True)
    run_simpleaf_quant(
        reads1=lib["reads1"],
        reads2=lib["reads2"],
        index=index,
        chemistry=chemistry,
        output=output,
        threads=threads,
        min_reads=min_reads,
        resolution=resolution,
        t2g_map=t2g,
        log_path=output.parent.parent / "logs" / f"quant_vdj_{receptor.lower()}.log",
    )
    h5 = find_quants_h5ad(output)
    adata = load_counts(h5, find_alevin_dir(output), f"vdj_{receptor}")
    ann = annotate_from_h5ad(adata, segments_tsv, receptor)
    csv_path = output / f"{receptor.lower()}_annotations.csv"
    ann.to_csv(csv_path, index=False)
    n_paired = int(ann["paired"].sum()) if "paired" in ann.columns else 0
    summary = {
        "receptor": receptor,
        "n_barcodes": int(adata.n_obs),
        "n_features": int(adata.n_vars),
        "n_paired": n_paired,
        "median_umi": float(np.median(ann["total_umi"])) if len(ann) else 0.0,
        "annotations": str(csv_path),
        "h5ad": str(h5),
    }
    (output / "vdj_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def run_vdj(cfg: dict, outdir: Path) -> dict | None:
    vdj = cfg.get("vdj")
    if not vdj:
        return None
    tcr = resolve_vdj_fastqs(vdj.get("tcr"))
    bcr = resolve_vdj_fastqs(vdj.get("bcr"))
    if not _reads_exist(tcr) and not _reads_exist(bcr):
        print("vdj: no TCR/BCR FASTQs found; skipping (GEX-only run)")
        return {"skipped": True, "reason": "FASTQs missing"}

    threads = int(vdj.get("threads", cfg.get("threads", 16)))
    chemistry = vdj.get("chemistry") or cfg.get("gex", {}).get("chemistry", "10xv3-5p")
    reference = vdj.get("reference")
    if not reference:
        raise ValueError("vdj.reference is required (10x refdata-cellranger-vdj directory or regions.fa)")
    idx = ensure_vdj_index(
        reference=Path(reference),
        index=Path(vdj["index"]) if vdj.get("index") else None,
        work_dir=outdir / "vdj",
        threads=threads,
        kmer_length=int(vdj.get("kmer_length", 21)),
        minimizer_length=int(vdj.get("minimizer_length", 11)),
    )
    t2g = idx / "t2g_3col.tsv"
    segments = Path(vdj.get("segments_tsv") or (idx.parent / "ref" / "segments.tsv"))
    if not segments.exists():
        raise FileNotFoundError(f"VDJ segments.tsv missing at {segments}")

    out = {"index": str(idx), "libraries": {}}
    if _reads_exist(tcr):
        out["libraries"]["TCR"] = run_vdj_library(
            lib=tcr,
            receptor="TCR",
            index=idx,
            t2g=t2g,
            segments_tsv=segments,
            chemistry=tcr.get("chemistry", chemistry),
            output=outdir / "vdj" / "tcr_quant",
            threads=int(tcr.get("threads", threads)),
            min_reads=int(tcr.get("min_reads", vdj.get("min_reads", 10))),
            resolution=tcr.get("resolution", vdj.get("resolution", "cr-like")),
        )
    if _reads_exist(bcr):
        out["libraries"]["BCR"] = run_vdj_library(
            lib=bcr,
            receptor="BCR",
            index=idx,
            t2g=t2g,
            segments_tsv=segments,
            chemistry=bcr.get("chemistry", chemistry),
            output=outdir / "vdj" / "bcr_quant",
            threads=int(bcr.get("threads", threads)),
            min_reads=int(bcr.get("min_reads", vdj.get("min_reads", 10))),
            resolution=bcr.get("resolution", vdj.get("resolution", "cr-like")),
        )
    (outdir / "vdj" / "vdj_run.json").write_text(json.dumps(out, indent=2) + "\n")
    return out


def attach_vdj_to_h5ad(h5ad_path: Path, vdj_run: dict | None) -> None:
    if not vdj_run or vdj_run.get("skipped") or not Path(h5ad_path).exists():
        return
    adata = ad.read_h5ad(h5ad_path)
    libs = vdj_run.get("libraries") or {}
    for receptor, spec in libs.items():
        csv_path = Path(spec["annotations"])
        if not csv_path.exists():
            continue
        ann = pd.read_csv(csv_path)
        prefix = "tcr" if receptor == "TCR" else "bcr"
        attach_vdj_obs(adata, ann, prefix)
    adata.write_h5ad(h5ad_path, compression="gzip")
