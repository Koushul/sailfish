# 48×48 microwell placement (entropy)

Python implementation of the localization used by
[https://lucid-crystal-kmqy.here.now/](https://lucid-crystal-kmqy.here.now/).

Each cell is placed on the chip by treating **row** and **column** independently.
Each axis has two plates of 48 spatial-hash oligos. Counts become a soft one-hot
multinomial log-likelihood (`β = 3`), plates are combined with averaged empirical
log-priors, softmax → posterior, then:

| quantity | formula |
|---|---|
| MAP well | `argmax π_row` × `argmax π_col` (indices 1..48) |
| confidence | `max(π_row) · max(π_col)` |
| discrete entropy | `H(π_row) + H(π_col)` bits |
| spatial entropy | separable 1-D Gaussian (`σ = 1.5` wells), then scaled so uniform still scores `2·log2(48)` |

`verify` recomputes MAP / posteriors / discrete entropy from the LLs embedded in the
site's `cells.js` and checks exact / tight numeric parity.

## Layout files

| file | role |
|---|---|
| `layouts/chip_layout.xlsx` | Source 384-well oligo plate (`Oligo Location via Plate`: 192 oligos, two plates) |
| `layouts/data.js` | Published `window.LAYOUT_DATA` used by the site (name, plate, axis, index, well, sequence) |

Plate 1 columns and plate 2 rows are reverse-numbered in the spreadsheet relative to
raw `oli010049–096` / `oli010145–192` order. **Always map ADT features through
`data.js` sequences**, not by `sbc*` name order.

ADT 15-mers live in `../refs/new_feature_ref_quant.csv` (`sbc1`–`sbc192`).

## Example

From the repo root (needs `numpy`; `from-h5ad` also needs `anndata`, `pandas`, `scipy`):

```bash
cd placement

python cell_placement.py selftest

# spreadsheet ↔ published data.js
python cell_placement.py verify-layout \
  --xlsx layouts/chip_layout.xlsx \
  --data-js layouts/data.js

# algorithm parity with the live E28S assignment on here.now
python cell_placement.py verify \
  --cells-js https://lucid-crystal-kmqy.here.now/datasets/E28S/cells.js

# place cells from an h5ad that has obsm['ADT'] (320 ADT columns in feature-ref order)
python cell_placement.py from-h5ad \
  --h5ad /path/to/sample_gex_adt.h5ad \
  --data-js layouts/data.js \
  --feature-ref ../refs/new_feature_ref_quant.csv \
  --min-layout-umi 10 \
  --out /tmp/assignments.csv
```

WagnerCollab MC38 (`mc38_velocity.h5ad`) has GEX only. Join Plate 1/2 ADT from
the E14S/E15S Cell Ranger matrices, localize with **Plate 1 only**, and write a
lucid-crystal-style site:

```bash
python export_wagner_plate1_site.py \
  --h5ad /ix1/ylee/shared/external/data/WagnerCollab/mc38_velocity.h5ad \
  --out sites/wagner_plate1
```

Rebuild `data.js` from the xlsx (should match the committed file):

```bash
python cell_placement.py layout-from-xlsx \
  --xlsx layouts/chip_layout.xlsx \
  --out /tmp/data.js
```
