import unittest
from pathlib import Path

from lib.vdj import (
    CHEMISTRY,
    GEOMETRY,
    TSO,
    annotate_read,
    call_cells,
    correct_barcode,
    find_cdr3,
    flank_score,
    parse_regions_fa,
    r1_insert,
    resolve_regions_fa,
)

EXAMPLE = (
    "AGACATCAGTGTACTTCTGTGCCAGCGGTGATGCGGGCAAAAACTCCGACTACACCTTCGGCTCAGGGACC"
    "AGGCTTTTGGTAATAGAGGATCTGAGAAATGTGACTCC"
)
REF = Path("/ix1/ylee/Yifan_Zhang/Cellranger_ref/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0")


class VdjTests(unittest.TestCase):
    def test_gemx_5p_geometry(self):
        self.assertEqual(GEOMETRY, "1{b[16]u[12]x:}2{r:}")
        self.assertEqual(CHEMISTRY, "10xv3-5p")
        r1 = "A" * 16 + "C" * 12 + TSO + EXAMPLE
        bc, umi, insert = r1_insert(r1)
        self.assertEqual(bc, "A" * 16)
        self.assertEqual(umi, "C" * 12)
        self.assertTrue(insert.startswith("AGACATCA"))

    def test_barcode_correction(self):
        wl = {"AACCGGTTAACCGGTT", "TTTTAAAACCCCGGGG"}
        self.assertEqual(correct_barcode("AACCGGTTAACCGGTT", wl), "AACCGGTTAACCGGTT")
        self.assertEqual(correct_barcode("AACCGGTTAACCGGTA", wl), "AACCGGTTAACCGGTT")
        self.assertIsNone(correct_barcode("AACCGGTTAACCGGCC", wl))

    def test_cdr3_fgxg(self):
        hit = find_cdr3(EXAMPLE)
        self.assertIsNotNone(hit)
        self.assertEqual(hit.aa, "CASGDAGKNSDYTF")
        self.assertTrue(hit.nt.startswith("TGT"))
        self.assertGreaterEqual(len(hit.nt), 26)
        self.assertLess(len(hit.nt), 80)
        self.assertGreaterEqual(hit.flank, 0)

    def test_flank_motifs_score_documented_example(self):
        prot = "LQPEDSAVYY" + "CASSL" + "LTFGAGTRVTV"
        score = flank_score(prot, 10, 15)
        self.assertGreaterEqual(score, 10)

    def test_short_v_overhang_still_annotates(self):
        if not REF.exists():
            self.skipTest("mouse V(D)J reference not on disk")
        from lib.vdj import translate

        ref = parse_regions_fa(resolve_regions_fa(REF))
        v = next(
            s
            for s in ref.v_by_chain["TRB"]
            if s.anchor and s.anchor >= 30 and translate(s.seq[s.anchor - 3 : s.anchor]) in "YFL"
        )
        j = next(s for s in ref.j_by_chain["TRB"] if s.anchor is not None and len(s.seq) >= s.anchor + 12)
        up = v.seq[v.anchor - 9 : v.anchor]
        seq = up + v.seq[v.anchor : v.anchor + 3] + ("GCT" * 8) + j.seq[j.anchor :]
        self.assertLess(len(up), 30)
        hit = find_cdr3(seq)
        self.assertIsNotNone(hit)
        self.assertLess(hit.start, 30)
        ann = annotate_read(seq, hit, ref, {"TRB"})
        self.assertIsNotNone(ann)
        self.assertEqual(ann["chain"], "TRB")
        self.assertTrue(ann["productive"])

    def test_one_umi_contigs_are_not_a_cell(self):
        rows = []
        for i in range(4):
            rows.append(
                {
                    "barcode": "AAAAAAAGTAAAAAAA",
                    "productive": True,
                    "chain": "TRA",
                    "v_gene": "TRAV1",
                    "cdr3_nt": f"TGT{'A' * i}",
                    "umis": 1,
                    "filtered_umis": 1,
                    "max_umi_reads": 2,
                    "umi_depths": [2],
                }
            )
        rows.append(
            {
                "barcode": "AAAAAAACAAAAAAAC",
                "productive": True,
                "chain": "TRB",
                "v_gene": "TRBV1",
                "cdr3_nt": "TGTC",
                "umis": 3,
                "filtered_umis": 3,
                "max_umi_reads": 4,
                "umi_depths": [4, 4, 2],
            }
        )
        call_cells(rows)
        self.assertFalse(rows[0]["is_cell"])
        self.assertTrue(rows[4]["is_cell"])

    def test_mouse_vdj_reference_parses(self):
        if not REF.exists():
            self.skipTest("mouse V(D)J reference not on disk")
        ref = parse_regions_fa(resolve_regions_fa(REF))
        self.assertGreater(len(ref.v_by_chain["TRB"]), 10)
        self.assertGreater(len(ref.j_by_chain["TRB"]), 5)
        self.assertTrue(all(seg.seq.startswith("ATG") for seg in ref.v_by_chain["TRB"]))
        self.assertTrue(all(seg.anchor is not None for seg in ref.j_by_chain["IGH"]))


if __name__ == "__main__":
    unittest.main()
