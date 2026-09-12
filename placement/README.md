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

Rebuild `data.js` from the xlsx (should match the committed file):

```bash
python cell_placement.py layout-from-xlsx \
  --xlsx layouts/chip_layout.xlsx \
  --out /tmp/data.js
```

## 2-barcode chip (`layout_2d.csv`)

When each well is one **row** oligo and one **column** oligo (no plate combining), use:

| file | role |
|---|---|
| `layouts/layout_2d.csv` | Chip map: Plate 1 `ROW 1–48` × Plate 2 `COLUMN 1–48` (other oligos in the CSV are unused) |
| `cellplacement_2barcodes.py` | Placement: `from-h5ad`, `from-counts`, `selftest` |
| `build_2barcode_site.py` | Writes a here.now-ready folder (`index.html`, `data.js`, `cells.js`, `assignments.csv`) |
| `run_2barcode_site.sh` | Wrapper around `build_2barcode_site.py` |
| `site_template_2barcodes/index.html` | Viewer template (dataset name injected at build time) |

Required inputs for the site pipeline: AnnData with `obsm['ADT']`, the layout CSV, and an ADT feature-ref CSV (`name`, `sequence`). Optional `obs['sample']` (or `sample_id` / `batch`) becomes sample checkboxes.

```bash
cd placement
python cellplacement_2barcodes.py selftest

python cellplacement_2barcodes.py from-h5ad \
  --h5ad /path/to/gex_adt.h5ad \
  --layout layouts/layout_2d.csv \
  --feature-ref ../refs/new_feature_ref_quant.csv \
  --min-layout-umi 10 \
  --out /tmp/assignments.csv

python build_2barcode_site.py \
  --h5ad /path/to/gex_adt.h5ad \
  --layout layouts/layout_2d.csv \
  --feature-ref ../refs/new_feature_ref_quant.csv \
  --dataset MyExp \
  --out /tmp/myexp_site \
  --min-layout-umi 10
# optional: --publish   (needs here.now publish.sh / $HERENOW_PUBLISH)
```
