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

## 5' V(D)J (TCR / BCR)

GEM-X Single Cell 5' v3 (CG000733): Read 1 is 16 bp barcode + 12 bp UMI + TSO `TTTCTTATATGGG`. Barcode inclusion list is `3M-5pgex-jan-2023` (Single Cell 5' v3). Chemistry name `10xv3-5p`, geometry `1{b[16]u[12]x:}2{r:}`. Recommended cycles are 28 / 10 / 10 / 90; extra Read 1 cycles are TSO plus the start of the insert.

A VDJ-only config omits `gex`. `vdj.tcr` / `vdj.bcr` take `fastqs` or `reads1`/`reads2`. Assembly follows the Cell Ranger V(D)J steps: whitelist correction (one mismatch), per-barcode CDR3s from the conserved cysteine through FGXG/WGXG (26–80 nt), reference-guided V/J contigs, productive = V start codon through the end of J, in frame, no stop. Read 2 often starts in J, so a CDR3 with only a short V overhang is kept when that overhang matches the V anchor. Each UMI votes for one CDR3. A cell needs a productive contig with at least two UMIs, three UMIs that have two reads, and a UMI deeper than 3% of the library N50. Cells need a productive contig, at least three UMIs with two read pairs, and a per-UMI depth of at least 3% of the N50. Clonotypes are shared productive CDR3 nucleotide sequences.

```json
{
  "sample": "E20T",
  "output": "runs/E20T",
  "vdj": {
    "reference": "/path/to/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0",
    "whitelist": "/path/to/3M-5pgex-jan-2023.txt.gz",
    "chemistry": "10xv3-5p",
    "tcr": {"fastqs": "/path/to/E20T_TCR"}
  }
}
```

**Outputs:** quant → `{output}/quant/{sample}_gex[_adt].h5ad` + MTX. OCM → `{output}/ocm/per_sample_outs/<id>/` and `{sample}_gex_adt_ocm.h5ad`. VDJ → `{output}/vdj/tcr` and/or `bcr` (`filtered_contig_annotations.csv`, `clonotypes.csv`, `airr_rearrangement.tsv`, `filtered_contig.fasta`). `cellranger.*` adds `compare_to_cellranger.json`.

GEM-X OCM uses a 2 bp overhang at barcode positions 7–8 (`GT` OB1, `CA` OB2, `TC` OB3, `AG` OB4). After VDJ calling, split the existing tables without another FASTQ pass:

```bash
python run.py configs/e20_e24/ocm_split.json --split-vdj-ocm
```

Each well is written to `{vdj_dir}/{tcr|bcr}/per_ocm/<sample_id>/` with clonotypes rebuilt inside that well. `filtered_contig_annotations.ocm.csv` keeps the pooled clonotype id in `pooled_clonotype_id`. If the config sets `h5ad`, the same tables are also written as one AnnData file with one row per productive contig and dataset, OCM well, chain, and clonotype fields in `obs`. `X` is that contig's UMI count.

E20–E24 GEX uses the same GEM-X 5' chemistry (`10xv3-5p`, piscem geometry `chromium_v3_5p`) and the mouse GRCm39 splici index. `run_e20_e24_gex.sh` quantifies each library and writes an OCM matrix per well under `{output}/ocm/per_sample_outs/<sample_id>/`. Configs: `configs/e20_e24/*_gex.json`.

## Microwell placement

`placement/cell_placement.py` is the entropy localizer used by
https://lucid-crystal-kmqy.here.now/ (MAP well + discrete/spatial entropy).
Layout: `placement/layouts/chip_layout.xlsx` (source oligos) and
`placement/layouts/data.js` (published 48×48 barcode map). See `placement/README.md`.

