# Sailfish

simpleaf / alevin-fry 10x GEX (± ADT), optional Cell Ranger–compatible OCM demux.

Needs `simpleaf`, `piscem`, `alevin-fry` on `PATH` (or `export SAILFISH_CONDA=/path/to/conda_env`). Python: `anndata`, `numpy`, `pandas`, `scipy`. Splici indexing also needs `roers` (`cargo install roers`). Sibling `conda_env` / `cargo_tools` are auto-detected. `ALEVIN_FRY_HOME` defaults to `{output}/af_home`.

```bash
python run.py my_sample.json
python run.py my_sample.json -o /tmp/out --skip-quant
sbatch run_sbatch.sh configs/e28s_ocm.json
```

`mode: "quant"` is pooled; `"ocm"` splits by the 2 bp overhang at barcode positions 7–8 (GEM-X: **GT=OB1, CA=OB2, TC=OB3, AG=OB4**).

## Indices

Config `*.index` must be the folder `simpleaf index --output DIR` writes as `DIR/index` (`piscem_idx*`, `t2g_3col.tsv`). Use a local `--work-dir`.

GEX splici (genome+GTF → `spliced` / `unspliced` / `ambiguous`). `--rlen` ≈ cDNA length (GEM-X v4 ~90; 10xv3 often 91). Default `--ref-type` is `spliced+intronic`.

```bash
export ALEVIN_FRY_HOME="${ALEVIN_FRY_HOME:-./af_home}"
simpleaf set-paths
simpleaf index \
  --output mouse-2024-A_splici \
  --fasta /path/to/refdata-gex-GRCm39-2024-A/fasta/genome.fa \
  --gtf /path/to/refdata-gex-GRCm39-2024-A/genes/genes.gtf \
  --ref-type spliced+intronic --rlen 91 --threads 16 \
  --work-dir /tmp/${USER}_simpleaf_index
# gex.index = mouse-2024-A_splici/index
```

ADT from a 10x feature CSV (`id`, `name`, `sequence`). Repo copy: `refs/new_feature_ref_quant.csv` (15 bp spatial + antibody). k=7, minimizer 5.

```bash
simpleaf index \
  --output adt_feature_index \
  --feature-csv refs/new_feature_ref_quant.csv \
  --kmer-length 7 --minimizer-length 5 --overwrite --threads 16
# adt.index = adt_feature_index/index
```

## Config

Omit `adt` for GEX-only; omit `ocm` (or `"mode": "quant"`) to skip demux. Default OCM map is OB1–OB4 if `ocm.samples` is missing.

```json
{
  "sample": "E28S",
  "mode": "ocm",
  "output": "runs/E28S",
  "threads": 16,
  "min_gex_umi": 500,
  "gex": {
    "reads1": ["fastqs/E28S_GEX_S3_R1_001.fastq.gz"],
    "reads2": ["fastqs/E28S_GEX_S3_R2_001.fastq.gz"],
    "index": "mouse-2024-A_splici/index",
    "chemistry": "10xv4-3p"
  },
  "adt": {
    "reads1": ["fastqs/E28S_ADT_S4_R1_001.fastq.gz"],
    "reads2": ["fastqs/E28S_ADT_S4_R2_001.fastq.gz"],
    "index": "adt_feature_index/index",
    "chemistry": "e14s-adt-10xv4",
    "geometry": "1{b[16]u[12]x:}2{r[15]x:}",
    "feature_ref": "refs/new_feature_ref_quant.csv"
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

| Field | Default | Meaning |
|-------|---------|---------|
| `sample` | `output` basename | Name in output files |
| `mode` | `quant` | `quant` pooled; `ocm` demux by overhang |
| `output` | required | Run dir |
| `threads` | `16` | simpleaf threads (GEX; ADT ~¼ unless `adt.threads`) |
| `skip_quant` | `false` | Skip mapping; reuse `gex.h5ad` or `{output}/gex_quant` |
| `min_gex_umi` | `500` | Cell filter (`filter.min_gex_umi` still works) |
| `alevin_fry_home` | `{output}/af_home` | Chemistry registry + permit lists |
| `tools.conda` | `$SAILFISH_CONDA` or sibling `conda_env` | Env with simpleaf/piscem/alevin-fry |
| `tools.cargo` | `$SAILFISH_CARGO` or sibling `cargo_tools` | Extra `bin/` (e.g. `roers`) |
| `gex.reads1` / `reads2` | | R1/R2 FASTQs (string or list, paired order) |
| `gex.index` | | `…/index` from `simpleaf index` |
| `gex.chemistry` | `10xv3` | e.g. `10xv3`, `10xv4-3p` |
| `gex.min_reads` | `10` | `--min-reads` with unfiltered permit list |
| `gex.resolution` | `cr-like` | UMI resolution (`cr-like`, `cr-like-em`, …) |
| `gex.h5ad` | | Existing GEX `quants.h5ad` (skip GEX map) |
| `gex.geometry` | | Register `gex.chemistry` with this geometry |
| `adt` | omitted | Same keys as `gex` |
| `adt.feature_ref` | | Feature CSV to rename ADT (`id`/`name`) |
| `adt.geometry` | | Custom chemistry (needed for `e14s-adt-10xv4`) |
| `ocm.samples` | GEM-X OB1–OB4 | Overhang → id `{"GT": "edge_gfp_plus"}`, or `{"OB1": "…"}`, or `[{sample_id, overhang}]` |
| `ocm.overhang_start` | `7` | 0-based overhang start in the 16 bp barcode |
| `ocm.overhang_len` | `2` | Overhang length |
| `ocm.include_unassigned` | `false` | Keep other overhangs as `unassigned` |
| `cellranger.per_sample_outs` | | `cellranger multi` `per_sample_outs` (OCM compare) |
| `cellranger.filtered_mtx` | | `cellranger count` `filtered_feature_bc_matrix` |
| `existing_quants.gex_h5ad` / `adt_h5ad` | | Same as `gex.h5ad` / `adt.h5ad` |

**Outputs:** quant → `{output}/quant/{sample}_gex[_adt].h5ad` + MTX. OCM → `{output}/ocm/per_sample_outs/<id>/` and `{sample}_gex_adt_ocm.h5ad`. `cellranger.*` adds `compare_to_cellranger.json`.
