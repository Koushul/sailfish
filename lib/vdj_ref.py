from __future__ import annotations

import json
from pathlib import Path

from .quant import run_simpleaf_index

REGION_ORDER = ("feature_id", "display_name", "gene_name", "region_type", "chain_type", "chain", "isotype", "allele")
KEEP_REGIONS = {"L-REGION+V-REGION", "V-REGION", "J-REGION", "C-REGION"}
DEFAULT_K = 21
DEFAULT_M = 11


def parse_vdj_header(line: str) -> dict:
    h = line[1:].strip() if line.startswith(">") else line.strip()
    parts = h.split("|")
    if len(parts) < 6:
        raise ValueError(f"unrecognized VDJ FASTA header: {line[:80]}")
    rec = {k: parts[i] if i < len(parts) else "" for i, k in enumerate(REGION_ORDER)}
    rec["feature_id"] = rec["feature_id"].split()[0]
    rec["gene_name"] = rec["gene_name"].strip()
    rec["locus"] = rec["chain"]
    rec["receptor"] = "TCR" if rec["chain_type"] == "TR" else "BCR" if rec["chain_type"] == "IG" else rec["chain_type"]
    return rec


def iter_regions_fa(path: Path):
    header = None
    seq: list[str] = []
    with Path(path).open() as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield parse_vdj_header(header), "".join(seq)
                header = line
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            yield parse_vdj_header(header), "".join(seq)


def resolve_regions_fa(reference: Path) -> Path:
    reference = Path(reference)
    if reference.is_file() and reference.name.endswith((".fa", ".fasta", ".fna")):
        return reference
    cand = reference / "fasta" / "regions.fa"
    if cand.exists():
        return cand
    raise FileNotFoundError(f"No fasta/regions.fa under {reference}")


def prepare_vdj_ref(reference: Path, outdir: Path, min_len: int = DEFAULT_K) -> dict:
    """Write simpleaf-ready FASTA + t2g from a 10x cellranger VDJ reference."""
    regions = resolve_regions_fa(reference)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fa_out = outdir / "vdj_segments.fa"
    t2g = outdir / "t2g_3col.tsv"
    meta = outdir / "segments.tsv"
    kept = []
    skipped_short = 0
    with fa_out.open("w") as fa, t2g.open("w") as tg, meta.open("w") as mt:
        mt.write("feature_id\tgene_name\tregion_type\tchain_type\tchain\tlocus\treceptor\tisotype\tallele\tlength\n")
        for rec, seq in iter_regions_fa(regions):
            if rec["region_type"] not in KEEP_REGIONS:
                continue
            if len(seq) < min_len:
                skipped_short += 1
                continue
            fid = rec["feature_id"]
            gene = rec["gene_name"].replace(" ", "_")
            fa.write(f">{fid}\n{seq}\n")
            tg.write(f"{fid}\t{gene}\tS\n")
            mt.write(
                f"{fid}\t{gene}\t{rec['region_type']}\t{rec['chain_type']}\t{rec['chain']}\t"
                f"{rec['locus']}\t{rec['receptor']}\t{rec['isotype']}\t{rec['allele']}\t{len(seq)}\n"
            )
            kept.append(rec)
    summary = {
        "source": str(regions),
        "n_segments": len(kept),
        "skipped_shorter_than_k": skipped_short,
        "min_len": min_len,
        "fasta": str(fa_out),
        "t2g": str(t2g),
        "segments_tsv": str(meta),
        "n_tcr": sum(1 for r in kept if rec_is_tcr(r)),
        "n_bcr": sum(1 for r in kept if rec_is_bcr(r)),
    }
    (outdir / "prepare_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if not kept:
        raise RuntimeError(f"No V/J/C segments ≥{min_len} bp in {regions}")
    return summary


def rec_is_tcr(rec: dict) -> bool:
    return rec.get("chain_type") == "TR" or str(rec.get("locus", "")).startswith("TR")


def rec_is_bcr(rec: dict) -> bool:
    return rec.get("chain_type") == "IG" or str(rec.get("locus", "")).startswith("IG")


def ensure_vdj_index(
    *,
    reference: Path,
    index: Path | None,
    work_dir: Path,
    threads: int,
    kmer_length: int = DEFAULT_K,
    minimizer_length: int = DEFAULT_M,
) -> Path:
    """Return path to simpleaf `…/index` directory, building it if needed."""
    if index is not None and Path(index).exists() and (Path(index) / "t2g_3col.tsv").exists():
        return Path(index)
    if index is not None and Path(index).exists() and list(Path(index).glob("piscem_idx*")):
        return Path(index)

    build_root = Path(index).parent if index is not None else Path(work_dir) / "vdj_index"
    prep = prepare_vdj_ref(Path(reference), build_root / "ref", min_len=kmer_length)
    run_simpleaf_index(
        ref_seq=Path(prep["fasta"]),
        output=build_root,
        threads=threads,
        kmer_length=kmer_length,
        minimizer_length=minimizer_length,
        keep_duplicates=True,
        work_dir=work_dir / "simpleaf_index_work",
    )
    idx = build_root / "index"
    t2g_idx = idx / "t2g_3col.tsv"
    t2g_idx.write_text(Path(prep["t2g"]).read_text())
    return idx
