# Sailfish

simpleaf / alevin-fry quantification for 10x GEX (± ADT), with optional Cell Ranger–compatible **OCM** demux.

## Setup

You need `simpleaf`, `piscem`, and `alevin-fry` on `PATH` (and a Python with `anndata`, `numpy`, `pandas`, `scipy`).

Either:

```bash
# conda env that already has simpleaf + piscem + alevin-fry
export SAILFISH_CONDA=/path/to/conda_env
```

or put those binaries on `PATH`. If this repo sits next to `af_tutorial/conda_env`, that env is picked up automatically.

```bash
python run.py my_sample.json
```

`ALEVIN_FRY_HOME` is created under `{output}/af_home` unless you set `alevin_fry_home`.

## Commands

```bash
# pooled GEX (± ADT)
python run.py configs/e28s_quant.json

# OCM demux (same quant, then split by barcode overhang)
python run.py configs/e28s_ocm.json

# override output dir; skip mapping if quants already exist
python run.py my_sample.json -o /tmp/out --skip-quant
```

Slurm: `sbatch run_sbatch.sh configs/e28s_ocm.json`

## Example config

Minimal **quant** (no OCM):

```json
{
  "sample": "pbmc",
  "output": "runs/pbmc",
  "gex": {
    "reads1": ["fastqs/sample_S1_L001_R1_001.fastq.gz"],
    "reads2": ["fastqs/sample_S1_L001_R2_001.fastq.gz"],
    "index": "/path/to/splici/index",
    "chemistry": "10xv3"
  }
}
```

Minimal **OCM** (GEM-X 4-plex defaults to GT/CA/TC/AG → OB1–OB4 unless `ocm.samples` is set):

```json
{
  "sample": "E28S",
  "mode": "ocm",
  "output": "runs/E28S",
  "gex": {
    "reads1": ["fastqs/E28S_GEX_S3_R1_001.fastq.gz"],
    "reads2": ["fastqs/E28S_GEX_S3_R2_001.fastq.gz"],
    "index": "/path/to/mouse-splici/index",
    "chemistry": "10xv4-3p"
  },
  "adt": {
    "reads1": ["fastqs/E28S_ADT_S4_R1_001.fastq.gz"],
    "reads2": ["fastqs/E28S_ADT_S4_R2_001.fastq.gz"],
    "index": "/path/to/adt_index",
    "chemistry": "e14s-adt-10xv4",
    "geometry": "1{b[16]u[12]x:}2{r[15]x:}",
    "feature_ref": "/path/to/feature_ref.csv"
  },
  "ocm": {
    "samples": {
      "GT": "edge_gfp_plus",
      "CA": "edge_gfp_minus",
      "TC": "core_gfp_plus",
      "AG": "core_gfp_minus"
    }
  }
}
```

### Parameters

| Field | Default | Meaning |
|-------|---------|---------|
| `sample` | basename of `output` | Name used in output files |
| `mode` | `"quant"` | `"quant"` = pooled cells. `"ocm"` = demux by 2 bp overhang in the 16 bp cell barcode |
| `output` | *(required)* | Run directory (quants, h5ad, logs) |
| `threads` | `16` | Threads for simpleaf (GEX uses this; ADT uses ~1/4 unless `adt.threads` is set) |
| `skip_quant` | `false` | If true, do not map; reuse `gex.h5ad` / previous `output/gex_quant` |
| `min_gex_umi` | `500` | Keep cells with at least this many GEX UMIs (`filter.min_gex_umi` still works) |
| `alevin_fry_home` | `{output}/af_home` | Chemistry registry + permit lists for simpleaf |
| `tools.conda` | `$SAILFISH_CONDA` or sibling `conda_env` | Conda prefix with simpleaf/piscem/alevin-fry |
| `tools.cargo` | `$SAILFISH_CARGO` or sibling `cargo_tools` | Extra `bin/` (e.g. `roers`) |
| **`gex.reads1` / `reads2`** | | R1/R2 FASTQs (string or list; order must match) |
| `gex.index` | | piscem/simpleaf index directory |
| `gex.chemistry` | `"10xv3"` | Registered chemistry (`10xv3`, `10xv4-3p`, …) |
| `gex.min_reads` | `10` | simpleaf `--min-reads` with unfiltered permit list |
| `gex.resolution` | `"cr-like"` | UMI resolution (`cr-like`, `cr-like-em`, …) |
| `gex.h5ad` | | Existing GEX `quants.h5ad` (implies skip mapping for GEX) |
| `gex.geometry` | | If set, register `gex.chemistry` with this geometry |
| **`adt`** | omitted | Same fields as `gex`. Omit the whole object for GEX-only |
| `adt.feature_ref` | | Cell Ranger–style feature CSV (id/name) to rename ADT features |
| `adt.geometry` | | Register a custom ADT chemistry (needed for `e14s-adt-10xv4`) |
| **`ocm.samples`** | GEM-X OB1–OB4 | Map overhang → sample id, e.g. `{"GT": "edge_gfp_plus", ...}`. Also accepts `{"OB1": "edge_gfp_plus"}` or a list of `{sample_id, overhang}` |
| `ocm.overhang_start` | `7` | 0-based start of the OCM overhang in the 16 bp barcode (Cell Ranger GEM-X) |
| `ocm.overhang_len` | `2` | Overhang length |
| `ocm.include_unassigned` | `false` | Keep barcodes whose overhang is not in `samples` as sample `unassigned` |
| `cellranger.per_sample_outs` | | Optional `cellranger multi` `outs/per_sample_outs` for OCM barcode comparison |
| `cellranger.filtered_mtx` | | Optional `cellranger count` `filtered_feature_bc_matrix` for pooled comparison |
| `existing_quants.gex_h5ad` / `adt_h5ad` | | Alternate to `gex.h5ad` / `adt.h5ad` |

GEM-X OCM overhangs: **GT=OB1, CA=OB2, TC=OB3, AG=OB4** (barcode bases 7–8).

## Outputs

- **quant:** `{output}/quant/{sample}_gex.h5ad` (or `_gex_adt.h5ad`) and 10x MTX folders
- **ocm:** `{output}/ocm/per_sample_outs/<sample_id>/` plus `{sample}_gex_adt_ocm.h5ad`

If `cellranger.*` is set, a `compare_to_cellranger.json` is written next to those outputs.
