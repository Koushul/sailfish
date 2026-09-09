#!/usr/bin/env python3
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from lib.vdj_ref import iter_regions_fa, parse_vdj_header, prepare_vdj_ref

MINI = ROOT / "tests" / "data" / "vdj_regions_mini.fa"
MOUSE_REF = Path(
    "/ix1/ylee/Yifan_Zhang/Cellranger_ref/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0"
)


class TestVdjRef(unittest.TestCase):
    def test_parse_header(self):
        rec = parse_vdj_header(">1|IGHA ENSMUST00000178282|IGHA|C-REGION|IG|IGH|A|00")
        self.assertEqual(rec["feature_id"], "1")
        self.assertEqual(rec["gene_name"], "IGHA")
        self.assertEqual(rec["region_type"], "C-REGION")
        self.assertEqual(rec["chain_type"], "IG")
        self.assertEqual(rec["locus"], "IGH")
        self.assertEqual(rec["receptor"], "BCR")

        rec = parse_vdj_header(">100|TRAV12-1 ENSMUST|TRAV12-1|L-REGION+V-REGION|TR|TRA|None|00")
        self.assertEqual(rec["receptor"], "TCR")
        self.assertEqual(rec["locus"], "TRA")

    def test_prepare_mini(self):
        with tempfile.TemporaryDirectory() as td:
            summary = prepare_vdj_ref(MINI, Path(td), min_len=21)
        self.assertEqual(summary["n_segments"], 3)
        self.assertEqual(summary["n_tcr"], 2)
        self.assertEqual(summary["n_bcr"], 1)
        self.assertGreaterEqual(summary["skipped_shorter_than_k"], 1)
        recs = list(iter_regions_fa(MINI))
        self.assertEqual(len(recs), 5)

    def test_annotate_dominant_vj(self):
        import anndata as ad
        import numpy as np
        from lib.vdj import annotate_from_h5ad

        with tempfile.TemporaryDirectory() as td:
            summary = prepare_vdj_ref(MINI, Path(td), min_len=21)
            adata = ad.AnnData(np.array([[3.0, 1.0, 0.0], [0.0, 0.0, 5.0]]))
            adata.obs_names = ["AAACCTGAAAAAAA-1", "AAACCTGACCCCCC-1"]
            adata.var_names = ["TRAV1", "TRAJ1", "IGHA"]
            tcr = annotate_from_h5ad(adata, Path(summary["segments_tsv"]), "TCR")
        self.assertTrue(bool(tcr.loc[0, "TRA_v"]) and bool(tcr.loc[0, "TRA_j"]))
        self.assertEqual(tcr.loc[0, "TRA_v"], "TRAV1")
        self.assertEqual(tcr.loc[0, "TRA_j"], "TRAJ1")
        self.assertFalse(bool(tcr.loc[0, "paired"]))
        self.assertEqual(tcr.loc[1, "TRA_v"], "")

    def test_skip_missing_fastqs(self):
        from lib.vdj import run_vdj

        with tempfile.TemporaryDirectory() as td:
            cfg = {
                "vdj": {
                    "reference": str(MOUSE_REF),
                    "tcr": {"fastqs": str(Path(td) / "missing_tcr")},
                },
                "threads": 1,
            }
            out = run_vdj(cfg, Path(td))
        self.assertTrue(out["skipped"])

    def test_mouse_regions_fa(self):
        if not (MOUSE_REF / "fasta" / "regions.fa").exists():
            self.skipTest("10x mouse VDJ reference not on this machine")
        recs = list(iter_regions_fa(MOUSE_REF / "fasta" / "regions.fa"))
        self.assertEqual(len(recs), 655)
        types = {r["region_type"] for r, _ in recs}
        self.assertIn("C-REGION", types)
        self.assertIn("D-REGION", types)
        self.assertTrue(any(r["chain"] == "TRA" for r, _ in recs))
        self.assertTrue(any(r["chain"] == "IGH" for r, _ in recs))
        with tempfile.TemporaryDirectory() as td:
            summary = prepare_vdj_ref(MOUSE_REF, Path(td), min_len=21)
        self.assertGreater(summary["n_tcr"], 50)
        self.assertGreater(summary["n_bcr"], 50)
        self.assertEqual(summary["n_segments"], summary["n_tcr"] + summary["n_bcr"])


if __name__ == "__main__":
    unittest.main()
