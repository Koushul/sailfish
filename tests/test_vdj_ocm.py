import csv
import tempfile
import unittest
from pathlib import Path

from lib.vdj_ocm import split_library

SAMPLES = [
    {"ocm_barcode_id": "OB1", "sample_id": "kidney_cd4", "overhang": "GT", "description": "kidney"},
    {"ocm_barcode_id": "OB2", "sample_id": "spleen_cd4", "overhang": "CA", "description": "spleen"},
]


def _row(barcode, chain, cdr3_nt, contig):
    return {
        "barcode": barcode,
        "is_cell": True,
        "contig_id": contig,
        "high_confidence": True,
        "length": 100,
        "chain": chain,
        "v_gene": "TRAV1",
        "d_gene": "",
        "j_gene": "TRAJ1",
        "c_gene": "TRAC",
        "full_length": True,
        "productive": True,
        "cdr3": "CAAAF",
        "cdr3_nt": cdr3_nt,
        "reads": 4,
        "umis": 2,
        "raw_clonotype_id": "clonotypeP",
        "raw_consensus_id": "clonotypeP_consensus_TRA",
    }


class VdjOcmTests(unittest.TestCase):
    def test_wells_do_not_share_clonotypes(self):
        gt = "AAAAAAAGTAAAAAAA"
        ca = "AAAAAAACAAAAAAAC"
        self.assertEqual(gt[7:9], "GT")
        self.assertEqual(ca[7:9], "CA")
        other = "AAAAAAATTAAAAAAA"
        rows = [
            _row(gt, "TRA", "TGTAAAAAA", "c1"),
            _row(gt, "TRB", "TGTBBBBBB", "c2"),
            _row(ca, "TRA", "TGTAAAAAA", "c3"),
            _row(other, "TRA", "TGTCCCCCC", "c4"),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            lib = Path(tmp)
            with (lib / "filtered_contig_annotations.csv").open("w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
                w.writeheader()
                w.writerows(rows)
            with (lib / "all_contig_annotations.csv").open("w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=["barcode", "is_cell"])
                w.writeheader()
                for row in rows:
                    w.writerow({"barcode": row["barcode"], "is_cell": True})
            with (lib / "airr_rearrangement.tsv").open("w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=["sequence_id", "sequence"], delimiter="\t")
                w.writeheader()
                for row in rows:
                    w.writerow({"sequence_id": row["contig_id"], "sequence": "ATGC"})
            metrics = split_library(lib, "TCR", SAMPLES)
            by_id = {m["sample_id"]: m for m in metrics}
            self.assertEqual(by_id["kidney_cd4"]["number_of_cells"], 1)
            self.assertEqual(by_id["kidney_cd4"]["number_of_paired_cells"], 1)
            self.assertEqual(by_id["spleen_cd4"]["number_of_cells"], 1)
            self.assertEqual(by_id["spleen_cd4"]["number_of_clonotypes"], 1)
            self.assertEqual(by_id["unassigned"]["number_of_cells"], 1)
            kidney = (lib / "per_ocm" / "kidney_cd4" / "clonotypes.csv").read_text()
            spleen = (lib / "per_ocm" / "spleen_cd4" / "clonotypes.csv").read_text()
            self.assertIn("TRA:TGTAAAAAA", kidney)
            self.assertNotIn("TRB", spleen)
            annotated = (lib / "filtered_contig_annotations.ocm.csv").read_text()
            self.assertIn("kidney_cd4", annotated)
            self.assertIn("clonotypeP", annotated)


if __name__ == "__main__":
    unittest.main()
