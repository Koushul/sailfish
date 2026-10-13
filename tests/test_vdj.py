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
    merge_variants,
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

    def test_plasma_cell_variants_merge_and_cell_is_kept(self):
        bc = "AAAAAAAAGAAAAAAA"
        top = "TGTGCAAGAGAGGAGTATTACTATGACTACTGG"

        def contig(chain, cdr3, umis):
            return {
                "barcode": bc,
                "productive": True,
                "chain": chain,
                "v_gene": "IGHV1-80" if chain == "IGH" else "IGKV1",
                "cdr3_nt": cdr3,
                "umis": umis,
                "reads": umis * 3,
                "filtered_umis": umis,
                "max_umi_reads": 5,
                "umi_depths": [3] * umis,
            }

        rows = [contig("IGH", top, 3000), contig("IGK", "TGTCAGCAGTATAGTAAGCTTCCGTATACG", 2500)]
        for i in range(len(top)):
            alt = "A" if top[i] != "A" else "C"
            rows.append(contig("IGH", top[:i] + alt + top[i + 1 :], 1))
        rows.append(contig("IGH", "TGTGCAAGGACCTACCGGTACTTCGATGTCTGG", 6))
        rows.append(contig("IGL", "TGTGCCCTATGGTACAGCAACCATTTCGTG", 4))
        merged = merge_variants(rows)
        igh = [r for r in merged if r["chain"] == "IGH"]
        self.assertEqual(len(igh), 2)
        self.assertEqual(max(r["umis"] for r in igh), 3000 + len(top))
        call_cells(merged)
        self.assertTrue(all(r["is_cell"] for r in merged))
        passed = [r for r in merged if r["passes_filter"]]
        self.assertEqual({r["chain"] for r in passed}, {"IGH", "IGK"})
        self.assertEqual(len(passed), 2)

    def test_gex_anchor_rejects_vdj_only_barcodes(self):
        def tcr(bc):
            return [
                {
                    "barcode": bc,
                    "productive": True,
                    "chain": chain,
                    "v_gene": f"{chain}V1",
                    "cdr3_nt": "TGTGCCAGCAGCTTT",
                    "umis": 3,
                    "reads": 9,
                    "filtered_umis": 3,
                    "max_umi_reads": 4,
                    "umi_depths": [4, 3, 2],
                }
                for chain in ("TRA", "TRB")
            ]

        rows = tcr("AAAAAAAGTAAAAAAA") + tcr("CCCCCCCGTCCCCCCC")
        stats = call_cells(rows, gex_cells={"AAAAAAAGTAAAAAAA"})
        self.assertTrue(rows[0]["is_cell"])
        self.assertFalse(rows[2]["is_cell"])
        self.assertEqual(stats["vdj_cells"], 2)
        self.assertEqual(stats["vdj_cells_without_gex"], 1)

    def test_three_strong_alleles_on_one_locus_is_a_doublet(self):
        rows = [
            {
                "barcode": "GGGGGGGTCGGGGGGG",
                "productive": True,
                "chain": "TRB",
                "v_gene": f"TRBV{i}",
                "cdr3_nt": f"TGT{'ACG'[i] * 12}TTT",
                "umis": 5,
                "reads": 15,
                "filtered_umis": 5,
                "max_umi_reads": 4,
                "umi_depths": [4] * 5,
            }
            for i in range(3)
        ]
        stats = call_cells(rows)
        self.assertFalse(rows[0]["is_cell"])
        self.assertTrue(rows[0]["multiplet"])
        self.assertEqual(stats["multiplets"], 1)

    def test_well_ambient_contig_is_removed_but_own_clone_member_kept(self):
        lung_top = "TGTGCAAGAGAGGAGTATTACTATGACTACTGG"

        def igh(bc, cdr3, umis):
            return {
                "barcode": bc,
                "productive": True,
                "chain": "IGH",
                "v_gene": "IGHV1-80",
                "cdr3_nt": cdr3,
                "umis": umis,
                "reads": umis * 3,
                "filtered_umis": umis,
                "max_umi_reads": 4,
                "umi_depths": [3] * umis,
            }

        rows = [igh(f"EMPTYAAAG{i:02d}AAAAA", lung_top, 3) for i in range(40)]
        b_cell = "BCELLAAAGTTTTTTT"
        own = "TGTGCAAGGACCTACCGGTACTTCGATGTCTGG"
        rows += [igh(b_cell, own, 12), igh(b_cell, lung_top, 3), igh(b_cell, "TGTACGAGGGGGTTAAGGGCTTACTGG", 1)]
        plasma = "PLASMAAAGTTTTTTT"
        rows += [igh(plasma, lung_top, 900), igh(plasma, own, 2)]
        call_cells(rows, gex_cells={b_cell, plasma})
        kept = {(r["barcode"], r["cdr3_nt"]) for r in rows if r.get("passes_filter")}
        self.assertIn((b_cell, own), kept)
        self.assertNotIn((b_cell, lung_top), kept)
        self.assertIn((plasma, lung_top), kept)
        self.assertNotIn((plasma, own), kept)

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
