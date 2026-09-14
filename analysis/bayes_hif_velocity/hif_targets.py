#!/usr/bin/env python3
"""Curated HIF-α transcriptional targets used as the hypoxia velocity panel.

These are downstream HIF-1α / HIF-2α genes (HRE-driven), not HIF-α itself.
Spliced levels track whether the HIF program is on; unspliced lag tracks whether
transcription is moving toward or away from hypoxia.
"""
from __future__ import annotations

from dataclasses import dataclass

# Direct HIF-α targets with repeated independent support (ChIP and/or HRE
# reporter + hypoxia induction + loss on reoxygenation / HIF knockdown).
# Hallmark HYPOXIA (Liberzon et al. 2015, MSigDB M5891) is a hypoxia-up set,
# not a TF-target set; the `in_hallmark` flag records membership only.
#
# Primary sources:
#   Mole et al. 2009 JBC (HIF-1α/HIF-2α ChIP)
#   Benita et al. 2009 NAR (core HIF-1 response across cell types)
#   Semenza reviews of glycolytic HIF-1 targets
#   Liberzon et al. 2015 Cell Systems (HALLMARK_HYPOXIA)
CORE_HIFA_TARGETS: tuple[tuple[str, str, str, str, bool], ...] = (
    # human, mouse, role, why it belongs, in HALLMARK_HYPOXIA
    ("SLC2A1", "Slc2a1", "glycolysis", "GLUT1; canonical HIF-1 glucose uptake", True),
    ("SLC2A3", "Slc2a3", "glycolysis", "GLUT3; HIF-1 glucose uptake", True),
    ("HK2", "Hk2", "glycolysis", "hexokinase 2; HIF-1 glycolytic enzyme", True),
    ("PFKL", "Pfkl", "glycolysis", "phosphofructokinase; HIF-1", True),
    ("PFKP", "Pfkp", "glycolysis", "platelet PFK; HIF-1", True),
    ("PFKFB3", "Pfkfb3", "glycolysis", "PFKFB3; HIF-1 glycolytic flux", True),
    ("ALDOA", "Aldoa", "glycolysis", "aldolase A; HIF-1", True),
    ("ALDOC", "Aldoc", "glycolysis", "aldolase C; HIF-1", True),
    ("TPI1", "Tpi1", "glycolysis", "triose phosphate isomerase; HIF-1", True),
    ("GAPDH", "Gapdh", "glycolysis", "GAPDH; HIF-1", True),
    ("PGK1", "Pgk1", "glycolysis", "PGK1; textbook HRE / HIF-1", True),
    ("ENO1", "Eno1", "glycolysis", "enolase 1; HIF-1", True),
    ("ENO2", "Eno2", "glycolysis", "enolase 2; HIF-1", True),
    ("PKM", "Pkm", "glycolysis", "PKM2; direct HIF-1 target", False),
    ("LDHA", "Ldha", "glycolysis", "LDHA; canonical HIF-1 lactate enzyme", True),
    ("PDK1", "Pdk1", "mitochondria", "PDK1; HIF-1 block of PDH / TCA entry", True),
    ("SLC16A3", "Slc16a3", "transport", "MCT4; HIF-1 lactate export", False),
    ("VEGFA", "Vegfa", "angiogenesis", "VEGFA; HIF-1 and HIF-2", True),
    ("ADM", "Adm", "peptide", "adrenomedullin; HIF peptide target", True),
    ("ANGPTL4", "Angptl4", "angiogenesis", "ANGPTL4; HIF angiogenic target", True),
    ("CA9", "Car9", "pH", "CA9; strongest HIF-1 reporter gene", False),
    ("BNIP3", "Bnip3", "mitophagy", "BNIP3; HIF-1 mitochondrial target", False),
    ("BNIP3L", "Bnip3l", "mitophagy", "BNIP3L/NIX; HIF mitophagy", True),
    ("NDRG1", "Ndrg1", "stress", "NDRG1; HIF-1 stress target", True),
    ("DDIT4", "Ddit4", "signaling", "REDD1; HIF-1 mTOR brake", True),
    ("P4HA1", "P4ha1", "ECM", "P4HA1; HIF collagen hydroxylase", True),
    ("P4HA2", "P4ha2", "ECM", "P4HA2; HIF collagen hydroxylase", True),
    ("LOX", "Lox", "ECM", "lysyl oxidase; HIF ECM target", True),
    ("EGLN3", "Egln3", "feedback", "PHD3; HIF-induced negative feedback", False),
    ("BHLHE40", "Bhlhe40", "transcription", "DEC1; HIF-1 transcriptional target", True),
    ("CXCR4", "Cxcr4", "migration", "CXCR4; HIF-1 chemokine receptor", True),
    ("SERPINE1", "Serpine1", "ECM", "PAI-1; HIF-1", True),
)

# Not used as velocity genes: oxygen-sensing machinery and the TF itself.
EXCLUDED_MACHINERY = (
    "HIF1A",
    "EPAS1",
    "ARNT",
    "VHL",
    "EGLN1",
    "EGLN2",
    "EPO",
)

HALLMARK_URL = "https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_HYPOXIA"


@dataclass(frozen=True)
class HifTarget:
    human: str
    mouse: str
    role: str
    note: str
    in_hallmark: bool

    @property
    def glycolysis(self) -> bool:
        return self.role == "glycolysis"


def core_targets() -> list[HifTarget]:
    return [HifTarget(*row) for row in CORE_HIFA_TARGETS]


def core_human_symbols() -> list[str]:
    return [g.human for g in core_targets()]


def core_mouse_symbols() -> list[str]:
    return [g.mouse for g in core_targets()]
