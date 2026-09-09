from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
import scipy.io

from .io_utils import load_barcode_set, strip_gem


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 1.0
    u = len(a | b)
    return len(a & b) / u if u else 0.0


def _load_mtx_dir(d: Path):
    barcodes = [strip_gem(x) for x in _lines(d / "barcodes.tsv.gz")]
    feats = []
    feat_ids = []
    feat_path = d / "features.tsv.gz"
    if not feat_path.exists():
        feat_path = d / "genes.tsv.gz"
    for line in _lines(feat_path):
        cols = line.split("\t")
        feat_ids.append(cols[0])
        feats.append(cols[1] if len(cols) > 1 else cols[0])
    mtx_path = d / "matrix.mtx.gz"
    if mtx_path.exists():
        opener = gzip.open
    else:
        mtx_path = d / "matrix.mtx"
        opener = open
    with opener(mtx_path, "rb") as fh:
        X = scipy.io.mmread(fh)
    X = sparse.csr_matrix(X)
    # 10x MTX is genes x cells
    if X.shape == (len(feat_ids), len(barcodes)):
        X = X.T.tocsr()
    elif X.shape != (len(barcodes), len(feat_ids)):
        raise ValueError(f"MTX shape {X.shape} vs {len(barcodes)} x {len(feat_ids)}")
    id_map = dict(zip(barcodes, range(len(barcodes))))
    gene_map = {gid: i for i, gid in enumerate(feat_ids)}
    gene_map.update({n: i for i, n in enumerate(feats) if n not in gene_map})
    return set(barcodes), gene_map, X, feat_ids, barcodes


def _lines(path: Path) -> list[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def compare_ocm_to_cellranger_multi(
    simpleaf_ocm_outdir: Path,
    cellranger_per_sample_outs: Path,
    out_json: Path,
    sample_ids: list[str] | None = None,
) -> dict:
    if sample_ids is None:
        sample_ids = sorted(
            p.name for p in (simpleaf_ocm_outdir / "per_sample_outs").iterdir() if p.is_dir()
        )
    rows = []
    for s in sample_ids:
        sf_raw = simpleaf_ocm_outdir / "per_sample_outs" / s / "sample_raw_feature_bc_matrix" / "barcodes.tsv.gz"
        sf_filt = simpleaf_ocm_outdir / "per_sample_outs" / s / "sample_filtered_feature_bc_matrix" / "barcodes.tsv.gz"
        sf_match = simpleaf_ocm_outdir / "per_sample_outs" / s / "sample_filtered_feature_bc_matrix_cr_match" / "barcodes.tsv.gz"
        cr_raw = cellranger_per_sample_outs / s / "sample_raw_feature_bc_matrix" / "barcodes.tsv.gz"
        cr_filt = cellranger_per_sample_outs / s / "sample_filtered_feature_bc_matrix" / "barcodes.tsv.gz"
        if not cr_filt.exists():
            continue
        A_raw, B_raw = load_barcode_set(sf_raw), load_barcode_set(cr_raw) if cr_raw.exists() else set()
        A_filt, B_filt = load_barcode_set(sf_filt), load_barcode_set(cr_filt)
        row = {
            "sample_id": s,
            "sf_raw": len(A_raw),
            "cr_raw": len(B_raw),
            "raw_intersection": len(A_raw & B_raw),
            "raw_jaccard": jaccard(A_raw, B_raw),
            "sf_filtered": len(A_filt),
            "cr_filtered": len(B_filt),
            "filt_intersection": len(A_filt & B_filt),
            "filt_recall_of_cr": len(A_filt & B_filt) / len(B_filt) if B_filt else 0.0,
            "filt_precision": len(A_filt & B_filt) / len(A_filt) if A_filt else 0.0,
            "filt_jaccard": jaccard(A_filt, B_filt),
            "cr_filt_in_sf_raw": len(B_filt & A_raw) / len(B_filt) if B_filt else 0.0,
        }
        if sf_match.exists():
            M = load_barcode_set(sf_match)
            row["cr_match_n"] = len(M)
            row["cr_match_exact"] = M == B_filt
        rows.append(row)
    payload = {
        "kind": "cellranger_multi",
        "samples": rows,
        "all_cr_filtered_recovered_in_sf_raw": bool(rows) and all(r["cr_filt_in_sf_raw"] == 1.0 for r in rows),
        "all_cr_filtered_recalled": bool(rows) and all(r["filt_recall_of_cr"] == 1.0 for r in rows),
    }
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2) + "\n")
    print(pd.DataFrame(rows).to_string(index=False) if rows else "no overlapping sample ids")
    print(f"Wrote {out_json}")
    return payload


def compare_to_cellranger_count(
    simpleaf_mtx_dir: Path,
    cellranger_mtx_dir: Path,
    out_json: Path,
    simpleaf_raw_barcodes: Path | None = None,
) -> dict:
    sf_bcs, sf_genes, sf_X, _, sf_bc_list = _load_mtx_dir(simpleaf_mtx_dir)
    cr_bcs, cr_genes, cr_X, cr_ids, cr_bc_list = _load_mtx_dir(cellranger_mtx_dir)
    shared_bc = sorted(sf_bcs & cr_bcs)
    shared_genes = sorted(set(sf_genes) & set(cr_genes))
    umi_corr = None
    gene_corr = None
    if shared_bc and shared_genes:
        sf_idx = {b: i for i, b in enumerate(sf_bc_list)}
        cr_idx = {b: i for i, b in enumerate(cr_bc_list)}
        sf_rows = [sf_idx[b] for b in shared_bc]
        cr_rows = [cr_idx[b] for b in shared_bc]
        gene_keys = shared_genes[: min(5000, len(shared_genes))]
        sf_cols = [sf_genes[g] for g in gene_keys]
        cr_cols = [cr_genes[g] for g in gene_keys]
        sf_sub = sf_X[sf_rows][:, sf_cols]
        cr_sub = cr_X[cr_rows][:, cr_cols]
        sf_umi = np.asarray(sf_sub.sum(1)).ravel()
        cr_umi = np.asarray(cr_sub.sum(1)).ravel()
        if sf_umi.std() > 0 and cr_umi.std() > 0:
            umi_corr = float(np.corrcoef(sf_umi, cr_umi)[0, 1])
        sf_g = np.asarray(sf_sub.sum(0)).ravel()
        cr_g = np.asarray(cr_sub.sum(0)).ravel()
        if sf_g.std() > 0 and cr_g.std() > 0:
            gene_corr = float(np.corrcoef(sf_g, cr_g)[0, 1])

    raw_recall = None
    if simpleaf_raw_barcodes is not None and simpleaf_raw_barcodes.exists():
        raw = load_barcode_set(simpleaf_raw_barcodes)
        raw_recall = len(cr_bcs & raw) / len(cr_bcs) if cr_bcs else 0.0

    payload = {
        "kind": "cellranger_count",
        "sf_filtered": len(sf_bcs),
        "cr_filtered": len(cr_bcs),
        "filt_intersection": len(sf_bcs & cr_bcs),
        "filt_recall_of_cr": len(sf_bcs & cr_bcs) / len(cr_bcs) if cr_bcs else 0.0,
        "filt_precision": len(sf_bcs & cr_bcs) / len(sf_bcs) if sf_bcs else 0.0,
        "filt_jaccard": jaccard(sf_bcs, cr_bcs),
        "n_shared_genes_for_corr": len(shared_genes),
        "cell_umi_pearson": umi_corr,
        "gene_umi_pearson": gene_corr,
        "cr_filt_in_sf_raw": raw_recall,
    }
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return payload
