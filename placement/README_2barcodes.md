# 2-barcode 48×48 site

Generate a static HTML viewer (same style as [lucid-crystal](https://lucid-crystal-kmqy.here.now/)) for chips that use **one row oligo × one column oligo** per well. No plate combining.

Python: `numpy`. Placement from h5ad and HTML export also need `anndata`, `pandas`, `scipy`.

## Required inputs

| input | what it is |
|---|---|
| `--h5ad` | AnnData with `obsm['ADT']`. Feature names from `uns['ADT_var']['feature_name']`, or pass `--adt-names`. Optional `obs['sample']` / `obs['sample_id']` / `obs['batch']` become sample checkboxes. |
| `--layout` | Chip CSV. Default: `layouts/layout_2d.csv`. |
| `--feature-ref` | ADT table with `name` and `sequence` (15-mer). Default: `../refs/new_feature_ref_quant.csv`. |

`layout_2d.csv` columns:

- `Name`
- `Sequence`
- `96 Source Plate #` (`Plate 1` / `Plate 2`)
- `Chip Location (Row/Column)` (`ROW n` / `COLUMN n`)
- `Well Position in 384 Source Plate`

Default families: **Plate 1 ROW 1–48** × **Plate 2 COLUMN 1–48**. Other oligos in the CSV are unused.

## 1. Check placement

From the repo root:

```bash
python placement/cellplacement_2barcodes.py selftest
```

Optional: assignments CSV only (no HTML):

```bash
python placement/cellplacement_2barcodes.py from-h5ad \
  --h5ad /path/to/gex_adt.h5ad \
  --layout placement/layouts/layout_2d.csv \
  --feature-ref refs/new_feature_ref_quant.csv \
  --min-layout-umi 10 \
  --out /tmp/assignments.csv
```

## 2. Generate HTML

```bash
python placement/build_2barcode_site.py \
  --h5ad /path/to/gex_adt.h5ad \
  --layout placement/layouts/layout_2d.csv \
  --feature-ref refs/new_feature_ref_quant.csv \
  --dataset MyExp \
  --out /tmp/myexp_site \
  --min-layout-umi 10
```

Same thing via the wrapper (must be run from anywhere; paths are resolved from the repo):

```bash
placement/run_2barcode_site.sh \
  --h5ad /path/to/gex_adt.h5ad \
  --layout placement/layouts/layout_2d.csv \
  --feature-ref refs/new_feature_ref_quant.csv \
  --dataset MyExp \
  --out /tmp/myexp_site \
  --min-layout-umi 10
```

`--min-layout-umi` default is **10** (cells below that UMI total on the 96 layout oligos are dropped).

Useful flags:

| flag | default | meaning |
|---|---|---|
| `--dataset` | `cells` | id in the UI and folder `datasets/<id>/` |
| `--dataset-label` | same as `--dataset` | dropdown label |
| `--title` | `<dataset> · 48×48 2-barcode localization` | browser title |
| `--heading` | `48×48 Microwell · 2-barcode layout` | page heading |
| `--subtitle` | auto (n cells, wells, UMI cutoff) | header subtitle (HTML allowed) |
| `--note` | `2 oligos/well · Plate-1 rows × Plate-2 columns` | legend note |
| `--row-plate` / `--row-axis` | `1` / `row` | which layout family is the row barcode |
| `--col-plate` / `--col-axis` | `2` / `column` | which layout family is the column barcode |
| `--adt-names` | (from h5ad) | text file, one ADT feature name per column |
| `--publish` | off | upload the output folder with here.now |

## 3. Output folder

```
/tmp/myexp_site/
  index.html                 # viewer (open this)
  data.js                    # layout oligos
  assignments.csv            # MAP well, entropy, confidence, layout UMI
  datasets/MyExp/cells.js    # per-cell posteriors + log-likelihoods
```

Preview locally:

```bash
python -m http.server 8000 --directory /tmp/myexp_site
# then open http://127.0.0.1:8000/
```

Do not open `index.html` as a `file://` URL if the browser blocks loading `cells.js`.

## 4. Publish to here.now (optional)

The builder looks for `publish.sh` at `$HERENOW_PUBLISH`, then `~/.claude/skills/here-now/scripts/publish.sh` or `~/.cursor/skills/here-now/scripts/publish.sh`.

```bash
python placement/build_2barcode_site.py \
  --h5ad /path/to/gex_adt.h5ad \
  --dataset MyExp \
  --out /tmp/myexp_site \
  --publish
```

Or publish a folder you already built:

```bash
"$HOME/.claude/skills/here-now/scripts/publish.sh" /tmp/myexp_site \
  --client cursor \
  --title "MyExp 48x48 2-barcode localization"
```

Authenticated publishes are permanent. The script prints the live URL (`https://….here.now/`).

## What the page shows

- **Barcodes** mode: 48×48 chip; click a well for the row oligo (Plate 1) and column oligo (Plate 2).
- **Cells** mode: MAP wells colored by entropy / count / confidence; click a cell for row/column posterior bars.
- Sample checkboxes appear only if the h5ad has `sample`, `sample_id`, or `batch`.
- **Spatial entropy** recomputes in the browser from the embedded row/column log-likelihoods.
