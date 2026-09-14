# A223 neutrophil HIF θ

Script: `python analysis/a223_neutrophil_hypoxia/run.py`

Finds neutrophils in sailfish E27/E29 OCM GEX (annotated Neutrophil labels + S100a8/Cxcr2 markers) and scores the E14/E15 Tumor HIF-down gene set. **A223 `hypoxia+` is Image-iT LIVE Green ROS (DCF), same FITC-band probe as E15**, used as a hypoxia surrogate; lactate+ is a second stain; DP = both. Persistent vs reverted gates use E14 neutrophil μ/σ (same 14 genes / weights as Tumor θ).

Results: `/ix1/ylee/kor11/A223/neutrophil_hypoxia/`

## Neutrophils in DCF+/hypoxia-sorted gates

There is **no GFP** in A223 or E14/E15 GEX. The FITC-band sort is Image-iT LIVE Green ROS (DCF). vs DN (not DCF+, not lactate+):

Annotated cells, Neutrophil labels:

| gate | n | % neutrophil | OR vs DN | p |
|------|--:|-------------:|---------:|--:|
| DN | 902 | **1.22%** | 1 | — |
| DCF+ (`hypoxia_plus`) | 14459 | **5.28%** | 4.52 | 4.9e-10 |
| DP | 12794 | 2.45% | 2.03 | 0.017 |
| DCF+ or DP | 27253 | 3.95% | 3.33 | 2.5e-6 |
| all annotated | 28410 | 3.85% | — | — |

DN matches unsorted A223 PMN-MDSC ~1.22% (Strait et al.). DCF+ over-represents neutrophils ~4.5-fold vs DN. E15S vs E14S placed cells: neutrophils **14.9% vs 8.3%** (OR 1.94).

Table: `neutrophil_gate_enrichment.csv`.

## Are DCF+ neutrophils making ROS according to RNA?

No. DCF is the ROS assay (oxidase assembly / p47phox phosphorylation). scRNA-seq does not capture that.

QC neutrophils (spliced UMI ≥ 1500): DCF+ n=272, DP n=111, DN n=4 (DN too small to test).

| module | DCF+ vs DP AUROC | interpretation |
|--------|-----------------:|----------------|
| NADPH oxidase subunits (`Cybb`, `Ncf1`, …) | 0.61 | small increase; `Cybb` still mostly undetected (28% vs 14% nonzero) |
| NRF2 / ROS response (`Hmox1`, `Nqo1`, `Sod2`, …) | 0.50 | no transcriptional ROS response |
| Wright priming (`Cxcl1/2`, `Il1b`, …) | 0.47 | not DCF-specific; DP is as high or higher |
| TNF/NF-κB priming | 0.41 | **higher in DP**, not DCF+ |
| HIF glycolytic module | 0.43 | **higher in DP** |

E15 vs E14 neutrophils (same genes): priming **up** (AUROC 0.69), oxidase **down** (0.40), NRF2 not up (0.46). Hypoxia-sorted neutrophils look primed, not oxidase-high.

Tables: `ros_rna_module_contrasts.csv`, `ros_rna_gene_medians.csv`.

## Dye vs true hypoxia (DCF+ neutrophils)

Hypothesis: neutrophils are FITC+/“hypoxia+” because they oxidize/handle Image-iT DCF, not because they sit in a HIF-high niche.

| test | dye prediction | result | dye? |
|------|----------------|--------|------|
| vs never-hypoxic E14 neu | DCF+ ≈ E14 | DCF+ 44% persist (θ 0.60) vs E14 21% (θ 0.95); AUROC 0.63 | partial — more HIF than E14, not fully hypoxic |
| DCF+-only vs DP (lactate) | DP is the hypoxic set | DP 62% persist (θ 0.04) vs DCF+ 44% (θ 0.60), p=1.3e-4 | **yes** |
| Same DCF+ gate, Tumor-scale θ | neus more reverted | Tumor 25% persist (θ 0.87); neu 24% persist (θ **0.97**, 67% reverted) | partial |
| Lineage vs HIF in DCF+ | neus enriched without HIF majority | neu OR 4.52 vs DN; only 44% of DCF+ neus HIF-persistent | **yes** |
| Mpo/esterase vs HIF | dye genes up, HIF not | Mpo/esterase AUROC 0.50; HIF AUROC 0.43 (higher in DP) | partial |

**Verdict:** DCF+ is a leaky hypoxia surrogate for neutrophils. Most DCF+-only neutrophils are HIF-reverted on the Tumor axis (median θ 0.97). HIF-persistent neutrophils concentrate in **DP** (DCF and lactate). Neutrophils are 4.5× over-represented in DCF+ anyway, which matches dye/ROS handling more than a HIF-high state. DCF+ is not identical to E14 (some extra HIF), so it is not pure false-positive dye with zero hypoxia biology.

`python analysis/a223_neutrophil_hypoxia/dye_vs_hypoxia.py` → `dye_vs_hypoxia_tests.csv`.
