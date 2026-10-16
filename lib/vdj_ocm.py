"""Split finished V(D)J tables by the GEM-X 2 bp overhang.

The overhang is already in the 16 bp barcode (0-based positions 7–8), so this
does not rescan FASTQs or recall cells. Clonotypes are rebuilt inside each
well from the productive contigs of cells that were already called.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

from .demux_ocm import overhang_of
from .vdj import CONTIG_COLS, _write_csv, assign_clonotypes

OCM_COLS = ["ocm_overhang", "ocm_barcode_id", "sample_id", *CONTIG_COLS, "pooled_clonotype_id"]


def _as_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "t", "yes"}


def _as_int(value) -> int:
    if value in (None, ""):
        return 0
    return int(float(value))


def _load_rows(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for key in ("is_cell", "high_confidence", "full_length", "productive"):
            if key in row:
                row[key] = _as_bool(row[key])
        for key in ("length", "reads", "umis"):
            if key in row:
                row[key] = _as_int(row[key])
        row.setdefault("sequence", "")
    return rows


def _load_sequences(libdir: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    airr = libdir / "airr_rearrangement.tsv"
    if airr.exists():
        with airr.open(newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                seqs[row["sequence_id"]] = row["sequence"]
        return seqs
    fasta = libdir / "filtered_contig.fasta"
    if not fasta.exists():
        return seqs
    name = ""
    chunks: list[str] = []
    with fasta.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name:
            seqs[name] = "".join(chunks)
    return seqs


def _assembled_barcodes(libdir: Path) -> dict[str, set[str]]:
    path = libdir / "all_contig_annotations.csv"
    if not path.exists():
        path = libdir / "filtered_contig_annotations.csv"
    found: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            bc = row["barcode"]
            found[overhang_of(bc)].add(bc)
    return found


def _paired(receptor: str, chains: set[str]) -> bool:
    if receptor == "TCR":
        return "TRA" in chains and "TRB" in chains
    if receptor == "BCR":
        return "IGH" in chains and bool({"IGK", "IGL"} & chains)
    return False


def _sample_index(samples: list[dict]) -> dict[str, dict]:
    return {s["overhang"]: s for s in samples}


def _prepare(row: dict, spec: dict | None) -> dict:
    out = dict(row)
    out["pooled_clonotype_id"] = row.get("raw_clonotype_id", "")
    out["raw_clonotype_id"] = ""
    out["raw_consensus_id"] = ""
    if spec is None:
        out["ocm_overhang"] = overhang_of(row["barcode"])
        out["ocm_barcode_id"] = "OTHER"
        out["sample_id"] = "unassigned"
    else:
        out["ocm_overhang"] = spec["overhang"]
        out["ocm_barcode_id"] = spec["ocm_barcode_id"]
        out["sample_id"] = spec["sample_id"]
    return out


def write_ocm_sample(
    outdir: Path,
    receptor: str,
    rows: list[dict],
    barcodes_assembled: int,
) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    rows.sort(key=lambda r: (r["barcode"], r["chain"], -r["umis"]))
    clonotypes, consensus = assign_clonotypes(rows)
    filtered = [r for r in rows if r.get("is_cell") and r.get("productive")]
    _write_csv(outdir / "filtered_contig_annotations.csv", filtered, OCM_COLS)
    _write_csv(
        outdir / "clonotypes.csv",
        clonotypes,
        ["clonotype_id", "frequency", "proportion", "cdr3s_aa", "cdr3s_nt"],
    )
    _write_csv(
        outdir / "consensus_annotations.csv",
        consensus,
        [
            "clonotype_id",
            "consensus_id",
            "length",
            "chain",
            "v_gene",
            "d_gene",
            "j_gene",
            "c_gene",
            "full_length",
            "productive",
            "cdr3",
            "cdr3_nt",
            "reads",
            "umis",
        ],
    )
    airr = []
    for row in filtered:
        airr.append(
            {
                "sequence_id": row["contig_id"],
                "sequence": row.get("sequence", ""),
                "productive": str(row["productive"]).upper(),
                "locus": row["chain"],
                "v_call": row["v_gene"],
                "d_call": row["d_gene"],
                "j_call": row["j_gene"],
                "c_call": row["c_gene"],
                "junction": row["cdr3_nt"],
                "junction_aa": row["cdr3"],
                "umi_count": row["umis"],
                "consensus_count": row["reads"],
                "cell_id": row["barcode"],
                "clone_id": row["raw_clonotype_id"],
            }
        )
    _write_csv(
        outdir / "airr_rearrangement.tsv",
        airr,
        [
            "sequence_id",
            "sequence",
            "productive",
            "locus",
            "v_call",
            "d_call",
            "j_call",
            "c_call",
            "junction",
            "junction_aa",
            "umi_count",
            "consensus_count",
            "cell_id",
            "clone_id",
        ],
        delimiter="\t",
    )
    with (outdir / "filtered_contig.fasta").open("w") as fh:
        for row in filtered:
            fh.write(f">{row['contig_id']}\n{row.get('sequence', '')}\n")
    by_bc: dict[str, set[str]] = defaultdict(set)
    for row in filtered:
        by_bc[row["barcode"]].add(row["chain"])
    paired = sum(1 for chains in by_bc.values() if _paired(receptor, chains))
    spec = rows[0] if rows else {}
    metrics = {
        "receptor": receptor,
        "sample_id": spec.get("sample_id", ""),
        "ocm_barcode_id": spec.get("ocm_barcode_id", ""),
        "ocm_overhang": spec.get("ocm_overhang", ""),
        "barcodes_assembled": barcodes_assembled,
        "number_of_cells": len(by_bc),
        "number_of_paired_cells": paired,
        "number_of_clonotypes": len(clonotypes),
        "number_of_productive_contigs": len(filtered),
    }
    with (outdir / "metrics_summary.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["metric", "value"])
        for key, value in metrics.items():
            w.writerow([key, value])
    return metrics


def split_library(libdir: Path, receptor: str, samples: list[dict]) -> list[dict]:
    filtered_path = libdir / "filtered_contig_annotations.csv"
    if not filtered_path.exists():
        raise FileNotFoundError(filtered_path)
    sequences = _load_sequences(libdir)
    rows = _load_rows(filtered_path)
    for row in rows:
        row["sequence"] = sequences.get(row["contig_id"], "")
    assembled = _assembled_barcodes(libdir)
    by_oh = _sample_index(samples)
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        oh = overhang_of(row["barcode"])
        spec = by_oh.get(oh)
        grouped[oh].append(_prepare(row, spec))

    per = libdir / "per_ocm"
    metrics = []
    annotated = []
    order = [s["overhang"] for s in samples]
    extras = [oh for oh in grouped if oh not in order]
    for oh in order + extras:
        group = grouped.get(oh, [])
        if not group and oh not in order:
            continue
        spec = by_oh.get(oh)
        sample_id = spec["sample_id"] if spec else "unassigned"
        if spec is None and not group:
            continue
        sample_rows = group
        if not sample_rows:
            sample_dir = per / sample_id
            sample_dir.mkdir(parents=True, exist_ok=True)
            _write_csv(sample_dir / "filtered_contig_annotations.csv", [], OCM_COLS)
            _write_csv(
                sample_dir / "clonotypes.csv",
                [],
                ["clonotype_id", "frequency", "proportion", "cdr3s_aa", "cdr3s_nt"],
            )
            hit = {
                "receptor": receptor,
                "sample_id": sample_id,
                "ocm_barcode_id": spec["ocm_barcode_id"] if spec else "OTHER",
                "ocm_overhang": oh if spec else "OTHER",
                "barcodes_assembled": len(assembled.get(oh, ())),
                "number_of_cells": 0,
                "number_of_paired_cells": 0,
                "number_of_clonotypes": 0,
                "number_of_productive_contigs": 0,
            }
            with (sample_dir / "metrics_summary.csv").open("w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["metric", "value"])
                for key, value in hit.items():
                    w.writerow([key, value])
            metrics.append(hit)
            continue
        hit = write_ocm_sample(per / sample_id, receptor, sample_rows, len(assembled.get(oh, ())))
        metrics.append(hit)
        annotated.extend(sample_rows)
    _write_csv(libdir / "filtered_contig_annotations.ocm.csv", annotated, OCM_COLS)
    return metrics


def split_vdj_ocm(cfg: dict) -> Path:
    libraries = cfg.get("libraries")
    if not libraries:
        raise ValueError("ocm split config needs 'libraries'")
    summary_rows = []
    out_summary = Path(cfg["summary"]) if cfg.get("summary") else None
    for lib in libraries:
        vdj_dir = Path(lib["vdj_dir"])
        samples = lib["samples"]
        (vdj_dir / "ocm_samples.json").write_text(json.dumps(samples, indent=2) + "\n")
        receptors = lib.get("receptors")
        if not receptors:
            receptors = [p.name for p in sorted(vdj_dir.iterdir()) if (p / "filtered_contig_annotations.csv").exists()]
        for receptor in receptors:
            libdir = vdj_dir / receptor
            print(f"ocm split {lib.get('id', vdj_dir.name)} {receptor}", flush=True)
            metrics = split_library(libdir, receptor.upper(), samples)
            for row in metrics:
                summary_rows.append(
                    {
                        "library": lib.get("id", vdj_dir.name),
                        "description": next(
                            (s.get("description", "") for s in samples if s["sample_id"] == row["sample_id"]),
                            "",
                        ),
                        **row,
                    }
                )
                print(
                    f"  {row['sample_id']} {row['ocm_overhang']}: "
                    f"{row['number_of_cells']} cells, {row['number_of_paired_cells']} paired, "
                    f"{row['number_of_clonotypes']} clonotypes",
                    flush=True,
                )
    fields = [
        "library",
        "receptor",
        "sample_id",
        "ocm_barcode_id",
        "ocm_overhang",
        "description",
        "barcodes_assembled",
        "number_of_cells",
        "number_of_paired_cells",
        "number_of_clonotypes",
        "number_of_productive_contigs",
    ]
    if out_summary is not None:
        out_summary.parent.mkdir(parents=True, exist_ok=True)
        _write_csv(out_summary, summary_rows, fields)
    result = out_summary or vdj_dir
    if cfg.get("h5ad"):
        result = write_vdj_h5ad(cfg, Path(cfg["h5ad"]))
    if cfg.get("joint_h5ad") or cfg.get("qc_cross_library"):
        from .vdj_joint import cross_library_qc, write_joint_h5ad

        if cfg.get("qc_cross_library"):
            cross_library_qc(cfg)
        if cfg.get("joint_h5ad"):
            write_joint_h5ad(cfg)
    return result


_TISSUES = ("submandibular", "lacrimal", "kidney", "spleen", "liver", "lung", "blood", "skin")


def tissue_population(sample_id: str) -> tuple[str, str]:
    for tissue in _TISSUES:
        prefix = tissue + "_"
        if sample_id.startswith(prefix):
            return tissue, sample_id[len(prefix) :]
    return "", sample_id


def _paired_cell(receptor: str, chains: set[str]) -> bool:
    if receptor == "TCR":
        return "TRA" in chains and "TRB" in chains
    if receptor == "BCR":
        return "IGH" in chains and bool({"IGK", "IGL"} & chains)
    return False


def write_vdj_h5ad(cfg: dict, path: Path) -> Path:
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse

    records = []
    sample_rows = []
    for lib in cfg["libraries"]:
        dataset = lib.get("id", Path(lib["vdj_dir"]).name)
        vdj_dir = Path(lib["vdj_dir"])
        by_id = {s["sample_id"]: s for s in lib["samples"]}
        for spec in lib["samples"]:
            tissue, population = tissue_population(spec["sample_id"])
            sample_rows.append(
                {
                    "dataset": dataset,
                    "sample_id": spec["sample_id"],
                    "description": spec.get("description", ""),
                    "ocm_barcode_id": spec["ocm_barcode_id"],
                    "ocm_overhang": spec["overhang"],
                    "tissue": tissue,
                    "population": population,
                    "sample_ul": float(spec["sample_ul"]) if "sample_ul" in spec else float("nan"),
                    "diluent_ul": float(spec["diluent_ul"]) if "diluent_ul" in spec else float("nan"),
                    "cells_sorted": float(spec["cells_sorted"]) if "cells_sorted" in spec else float("nan"),
                }
            )
        for receptor_name in lib["receptors"]:
            receptor = receptor_name.upper()
            libdir = vdj_dir / receptor_name
            for spec in lib["samples"]:
                sample_dir = libdir / "per_ocm" / spec["sample_id"]
                contigs = _load_rows(sample_dir / "filtered_contig_annotations.csv")
                sequences = _load_sequences(sample_dir)
                clones = {
                    row["clonotype_id"]: row
                    for row in csv.DictReader((sample_dir / "clonotypes.csv").open(newline=""))
                }
                tissue, population = tissue_population(spec["sample_id"])
                for row in contigs:
                    cid = row.get("raw_clonotype_id", "")
                    clone = clones.get(cid, {})
                    records.append(
                        {
                            "dataset": dataset,
                            "receptor": receptor,
                            "organism": "Mus musculus",
                            "chemistry": "10xv3-5p",
                            "barcode": row["barcode"],
                            "cell_id": f"{dataset}|{receptor}|{row['barcode']}",
                            "contig_id": row["contig_id"],
                            "ocm_overhang": row["ocm_overhang"],
                            "ocm_barcode_id": row["ocm_barcode_id"],
                            "sample_id": row["sample_id"],
                            "description": spec.get("description", ""),
                            "tissue": tissue,
                            "population": population,
                            "sample_ul": float(spec["sample_ul"]) if "sample_ul" in spec else float("nan"),
                            "diluent_ul": float(spec["diluent_ul"]) if "diluent_ul" in spec else float("nan"),
                            "cells_sorted": float(spec["cells_sorted"]) if "cells_sorted" in spec else float("nan"),
                            "chain": row["chain"],
                            "v_gene": row["v_gene"],
                            "d_gene": row["d_gene"],
                            "j_gene": row["j_gene"],
                            "c_gene": row["c_gene"],
                            "cdr3": row["cdr3"],
                            "cdr3_nt": row["cdr3_nt"],
                            "sequence": sequences.get(row["contig_id"], row.get("sequence", "")),
                            "length": row["length"],
                            "reads": row["reads"],
                            "umis": row["umis"],
                            "full_length": row["full_length"],
                            "productive": row["productive"],
                            "high_confidence": row["high_confidence"],
                            "is_cell": row["is_cell"],
                            "well_clonotype_id": cid,
                            "clonotype_id": f"{dataset}|{receptor}|{spec['sample_id']}|{cid}",
                            "pooled_clonotype_id": row.get("pooled_clonotype_id", ""),
                            "consensus_id": row.get("raw_consensus_id", ""),
                            "clonotype_frequency": int(clone.get("frequency", 0) or 0),
                            "clonotype_proportion": float(clone.get("proportion", 0) or 0),
                        }
                    )
            if by_id and not (libdir / "per_ocm").exists():
                raise FileNotFoundError(libdir / "per_ocm")
    if not records:
        raise ValueError("no VDJ contigs to write")
    obs = pd.DataFrame.from_records(records)
    chain_sets = obs.groupby("cell_id", sort=False)["chain"].agg(lambda s: set(s))
    receptor_of = obs.groupby("cell_id", sort=False)["receptor"].first()
    paired = {
        cell: _paired_cell(receptor_of.loc[cell], chains) for cell, chains in chain_sets.items()
    }
    obs["is_paired"] = obs["cell_id"].map(paired)
    obs["n_contigs"] = obs.groupby("cell_id", sort=False)["contig_id"].transform("size").astype(int)
    obs.index = pd.Index(
        obs["dataset"] + "|" + obs["receptor"] + "|" + obs["contig_id"],
        name="contig_key",
    )
    if not obs.index.is_unique:
        raise ValueError("contig keys are not unique")
    for col in (
        "dataset",
        "receptor",
        "organism",
        "chemistry",
        "ocm_overhang",
        "ocm_barcode_id",
        "sample_id",
        "description",
        "tissue",
        "population",
        "chain",
        "v_gene",
        "d_gene",
        "j_gene",
        "c_gene",
        "well_clonotype_id",
    ):
        obs[col] = obs[col].astype("category")
    x = sparse.csr_matrix(obs["umis"].to_numpy(np.float32).reshape(-1, 1))
    var = pd.DataFrame(index=pd.Index(["umis"], name="feature"))
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.uns["organism"] = "Mus musculus"
    adata.uns["chemistry"] = "10xv3-5p"
    adata.uns["barcode_whitelist"] = "3M-5pgex-jan-2023"
    adata.uns["ocm_overhang_start"] = 7
    adata.uns["ocm_overhang_len"] = 2
    adata.uns["ocm_overhangs"] = {"OB1": "GT", "OB2": "CA", "OB3": "TC", "OB4": "AG"}
    adata.uns["ocm_samples"] = pd.DataFrame(sample_rows)
    adata.uns["x"] = "UMI count of this productive contig"
    if cfg.get("summary") and Path(cfg["summary"]).exists():
        adata.uns["ocm_summary"] = pd.read_csv(cfg["summary"])
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path, compression="gzip")
    print(f"wrote {path} ({adata.n_obs} contigs, {obs['cell_id'].nunique()} cells)", flush=True)
    return path
