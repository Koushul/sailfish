# Assessment of HIF-α lag fits

What the fits actually support, ranked by how much the design can identify. Forced control persist is always 0% and is ignored below; empirical θ-gates are the measurement.

## Verdict

1. **The lag / velocity half of the model is not identified on any real library.** \(\sigma_v\) sits on the \(e^{0.5}\) clip in every fit. Hard flux is ~0. Tumor \(p^{\mathrm{away}}\) is never higher in the exposed arm than in the never-hypoxic (or DN) arm. There is no evidence of a reoxygenation or induction wave in unspliced HIF-target RNA.
2. **Spliced persist fractions are not a stable biological number.** They move with library-size ratio, cell typing, chemistry, and how small the control is. Do not quote a single “% persistent hypoxia.”
3. **The one hypoxia-exposure experiment that has a real never-hypoxic arm (E14 vs E15) does not show a large persistent HIF-on block** after CPM. On the placed object, E15 tumor persist (7.6%) is **below** the E14 empirical tail (11.2%). Wagner’s processing of the same libraries gives 20% — so even that experiment is processing-dependent.
4. **A223 uses the same Image-iT GFP-channel probe as E15S**, so dye+ vs HIF-target RNA is on-target. It is still not transgenic GFP and not HIF-α protein. E27 vs E29 persist does not replicate. DN n is 12–13 tumors. Neutrophil persist is a 3–4 cell MAD.

## 1. E14 vs E15 (only real never-hypoxic vs hypoxia-exposed design)

n_control is adequate for tumors (152) and neutrophils (319).

| | empirical control persist | exposed persist | UMI ratio | \(p^{\mathrm{away}}\) ctrl vs exp |
|--|---------------------------|-----------------|-----------|-------------------------------------|
| Tumor (placed) | 11.2% | **7.6%** | 4.77 | 0.16 vs 0.14 |
| Tumor (Wagner, same libraries) | 11.5% | **20.4%** | 1.01 | 0.083 vs 0.020 |
| neutrophil (placed) | 5.3% | 6.1% | 0.84 | 0.021 vs 0.011 |
| neutrophil (Wagner) | 4.6% | 3.4% | 0.98 | 0.061 vs 0.042 |

**Read this as:** a typical never-hypoxic tumor already has an ~11% tail above \(\theta=0.3\). That tail is reproducible (placed 11.2%, Wagner 11.5%, synthetics ~14%). E15 persist on the placed object does not exceed it. The historical 30% (placed) and 91% (Wagner `hypoxia_state`) persist calls mixed depth or an old score into the program; agreement with the new gate is 58% and 23%.

Wagner vs placed disagree on **exposed** persist because depth is matched in Wagner (median ~8.7k both samples) and 5× different in the placed object (14k vs 68k). CPM + within-sample \(\ell_n\) removes the 5× batch but also shrinks the E15 shift. Persist after hypoxia is **not identified across processings** of the same experiment.

Gene loadings on placed tumors are glycolytic (`Pgk1`, `Aldoa`, `Tpi1`, `Slc2a1`) more than CA9 (`Car9` \(\beta=0.29\) vs `Pgk1` \(0.49\)). Wagner tumors put more weight on `Slc2a1` / `Ndrg1` / `Bnip3`; `Car9` is still modest. That is a HIF-adjacent metabolic program, not a CA9-high hypoxia call.

Neutrophil *Car9* is off. The neutrophil factor is glycolysis (`Pkm`, `Aldoa`, `Eno1`). GMM collapsed to 100% reverted. Neutrophil persist ~3–6% is the θ-gate only, and \(h\) still tracks \(\log L\) on the placed object (Spearman −0.30).

## 2. A223 E27 / E29 (same Image-iT probe as E15S, OCM-gated)

`hypoxia_plus` is Image-iT+ in the GFP channel, **the same probe as E15S**. DN is Image-iT−. DP is Image-iT+ and lactate+. That is dye vs HIF-target RNA, not a ROS assay. Design still differs from E14/E15: E15S is a whole hypoxia-exposed library; A223 is OCM gates inside each chemistry.

Control is DN, n=13/12 tumors and 4/3 neutrophils. That is not a control MAD. Mean DN \(h\) on E29 neutrophils is −2.3 with exposed \(h\) clipped near 7 — the scale is exploding.

| tumor gate | E27 persist | E27 mean \(h\) | E29 persist | E29 mean \(h\) |
|------------|-------------|----------------|-------------|----------------|
| DN (Image-iT−) | 7.7% (n=13) | 0.12 | 0% (n=12) | 0.00 |
| hypoxia_plus (Image-iT+) | 35% | 1.03 | 6% | 0.12 |
| DP (Image-iT+ and lactate+) | 43% | 1.38 | 8% | 0.34 |

Within E27, mean \(h\) orders DP > Image-iT+ > DN, which is the expected dye/lactate stacking if the HIF-target factor tracks the same oxygen reporter as E15S. It does **not** replicate on E29 (both exposed gates stay near the DN location). UMI ratio is inverted across chemistry (E27 0.52, E29 2.75). *Car9* loads on E27 tumors (\(\beta=0.68\)) but is undetected for velocity on E29.

So: **E27 is compatible with Image-iT+ tumors having a higher HIF-target program than Image-iT−**; **E29 is not a confirmation**. Persist % itself remains unidentified because n_DN is 12–13 and 3′ vs 5′ disagree.

Neutrophils: 72–95% persist on Image-iT+/DP. That is the 3–4 cell control, plus neutrophil \(h\) correlating with \(\log L\) (+0.36 to +0.46). Ignore as a HIF-target fraction.

Tumor \(p^{\mathrm{away}}\) is ~0.01–0.02 on DN and Image-iT+. No lag difference by dye gate.

Do not call Image-iT+ transgenic GFP or HIF-α protein. Do not pool E27 with E29.

## 3. What is consistent across every fit

- Cycle is off \(\xi\) (Spearman \(\xi\) vs S near 0). The unspliced cycle covariate is doing its job; it is not why lag is empty.
- \(\sigma_v = e^{0.5}\) always. The lag scale prior/clip is binding. Unspliced HIF-target counts are too sparse/noisy for \(\xi\).
- Hard `transitioning_out` is not a usable label.
- Two phenotypes disagree (θ-gate vs GMM), worst in neutrophils.
- Balanced WLS still gives the tiny control 50% of the factor weight (A223 tumors: 13 cells vs ~5000).

## 4. What to report (and what not to)

**Report**

- Empirical never-hypoxic tail ~10% on MC38 tumors when n_control is tens to hundreds.
- After depth correction, E15 is not a large persistent HIF-on compartment on the placed object; historical persist % was confounded.
- Wagner shows the same tail and a larger E15 shift when depth is matched — processing sensitivity, not a second experiment.
- Lag does not add a direction of travel on these libraries.
- A223: same Image-iT probe as E15S; E27 Image-iT+ tumors have higher \(h\) than DN; E29 does not replicate; DN n too small for a persist %; neutrophil persist unusable.

**Do not report**

- A single persist percentage as the hypoxia fraction.
- E14 0% persist (forced).
- E15 reoxygenation from \(p^{\mathrm{away}}\).
- Neutrophil persist as shared tumor HIF kinetics.
- A223 Image-iT+ as transgenic GFP or as HIF-α protein; do not treat E27 persist % as confirmed by E29.
- GMM and θ-gate as independent confirmation.

The model is a **spliced HIF-target factor with a failed lag add-on**. Phenotype can be used, with the empirical control tail as the noise floor. Direction of hypoxia (toward/away) should not be claimed from these fits.
