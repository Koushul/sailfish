#!/bin/bash
# Quantify E20–E24 GEM-X 5' GEX and split each library by OCM well.
# Reads and the splici index are staged on node-local scratch. Final matrices
# are copied back next to the VDJ results. This does not rerun VDJ.

set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
PY="${PY:-python3}"
STAGE="${STAGE:-/scratch/slurm-4343839/e20_e24_gex}"
INDEX_SRC="/ix1/ylee/shared/E18S_E19S_E37S1_E37S2/E37S2_sailfish/refs/mouse-2024-A_splici/index"
RAW="/ix1/ylee/shared/E20_E30_Series/01.RawData"
DEST="/ix1/ylee/kor11/TCR_BCR/runs"
LOGDIR="/ix1/ylee/kor11/TCR_BCR/logs"
INDEX="$STAGE/index"
AFHOME="$STAGE/af_home"
mkdir -p "$STAGE/fastq" "$INDEX" "$LOGDIR"

if [[ ! -s "$INDEX/t2g_3col.tsv" ]]; then
  echo "STAGE INDEX $(date -Is)"
  rsync -a "$INDEX_SRC/" "$INDEX/"
fi

stage_pair() {
  local src_r1="$1" src_r2="$2" name="$3"
  local dst_r1="$STAGE/fastq/${name}_R1.fastq.gz"
  local dst_r2="$STAGE/fastq/${name}_R2.fastq.gz"
  if [[ ! -s "$dst_r1" ]]; then
    echo "COPY R1 $name $(date -Is)"
    cp -a "$src_r1" "$dst_r1"
  fi
  if [[ ! -s "$dst_r2" ]]; then
    echo "COPY R2 $name $(date -Is)"
    cp -a "$src_r2" "$dst_r2"
  fi
}

write_staged_config() {
  local cfg="$1" name="$2" out="$3"
  "$PY" - "$cfg" "$name" "$out" "$INDEX" "$AFHOME" "$STAGE" << 'PY'
import json, sys
from pathlib import Path
src, name, out, index, afhome, stage = sys.argv[1:]
cfg = json.loads(Path(src).read_text())
cfg["output"] = out
cfg["alevin_fry_home"] = afhome
cfg["gex"]["index"] = index
cfg["gex"]["reads1"] = [f"{stage}/fastq/{name}_R1.fastq.gz"]
cfg["gex"]["reads2"] = [f"{stage}/fastq/{name}_R2.fastq.gz"]
Path(out).mkdir(parents=True, exist_ok=True)
dest = Path(out) / "staged_config.json"
dest.write_text(json.dumps(cfg, indent=2) + "\n")
print(dest)
PY
}

copy_back() {
  local name="$1" src="$2"
  local dst="$DEST/$name"
  mkdir -p "$dst/gex_quant/af_quant/alevin" "$dst/logs"
  rsync -a "$src/gex_quant/af_quant/alevin/" "$dst/gex_quant/af_quant/alevin/"
  if [[ -d "$src/ocm" ]]; then
    rsync -a "$src/ocm/" "$dst/ocm/"
  fi
  if [[ -d "$src/logs" ]]; then
    rsync -a "$src/logs/" "$dst/logs/"
  fi
  cp -a "$src/run_config.json" "$dst/run_config_gex.json"
  if [[ ! -s "$dst/gex_quant/af_quant/alevin/quants.h5ad" ]]; then
    echo "missing quants.h5ad for $name" >&2
    exit 1
  fi
}

run_one() {
  local cfg="$1" name="$2" r1="$3" r2="$4"
  local out="$STAGE/out/$name"
  echo "START $(date -Is) $name"
  stage_pair "$r1" "$r2" "$name"
  local staged
  staged="$(write_staged_config "$cfg" "$name" "$out")"
  "$PY" "$ROOT/run.py" "$staged" > "$LOGDIR/${name}_gex.log" 2>&1
  copy_back "$name" "$out"
  rm -f "$STAGE/fastq/${name}_R1.fastq.gz" "$STAGE/fastq/${name}_R2.fastq.gz"
  if [[ -d "$out/gex_quant/af_map" ]]; then
    rm -rf "$out/gex_quant/af_map"
  fi
  echo "DONE $(date -Is) $name"
}

echo "PILOT $(date -Is)"
python3 - << PY
import gzip
from pathlib import Path
src_r1 = Path("$RAW/E20T_GEX/E20T_GEX-SI_TT_E4_23W5C2LT4_S1_L005_R1_001.fastq.gz")
src_r2 = Path("$RAW/E20T_GEX/E20T_GEX-SI_TT_E4_23W5C2LT4_S1_L005_R2_001.fastq.gz")
out = Path("$STAGE/fastq")
n = 40000
for src, dest in ((src_r1, out / "pilot_R1.fastq.gz"), (src_r2, out / "pilot_R2.fastq.gz")):
    with gzip.open(src, "rt") as inp, gzip.open(dest, "wt") as outp:
        for i, line in enumerate(inp):
            if i >= n * 4:
                break
            outp.write(line)
    print("wrote", dest, flush=True)
PY
pilot_out="$STAGE/out/pilot"
staged="$(write_staged_config "$ROOT/configs/e20_e24/e20t_gex.json" pilot "$pilot_out")"
"$PY" "$ROOT/run.py" "$staged" > "$LOGDIR/pilot_gex.log" 2>&1
python3 - << PY
import json
from pathlib import Path
info = json.loads(Path("$pilot_out/gex_quant/af_map/map_info.json").read_text())
pct = float(info["percent_mapped"])
print(f"pilot mapped {info['num_mapped']}/{info['num_reads']} ({pct}%) geometry check in cmdline", flush=True)
cmd = info.get("cmdline", "")
if "chromium_v3_5p" not in cmd:
    raise SystemExit(f"expected chromium_v3_5p in piscem command: {cmd}")
if pct < 50:
    raise SystemExit(f"mapping rate {pct} is below 50")
PY
rm -f "$STAGE/fastq/pilot_R1.fastq.gz" "$STAGE/fastq/pilot_R2.fastq.gz"
rm -rf "$pilot_out/gex_quant/af_map"
echo "PILOT OK $(date -Is)"

run_one "$ROOT/configs/e20_e24/e20t_gex.json" E20T \
  "$RAW/E20T_GEX/E20T_GEX-SI_TT_E4_23W5C2LT4_S1_L005_R1_001.fastq.gz" \
  "$RAW/E20T_GEX/E20T_GEX-SI_TT_E4_23W5C2LT4_S1_L005_R2_001.fastq.gz"
run_one "$ROOT/configs/e20_e24/e21t_gex.json" E21T \
  "$RAW/E21T_GEX/E21T_GEX-SI_TT_E5_23W5C2LT4_S2_L005_R1_001.fastq.gz" \
  "$RAW/E21T_GEX/E21T_GEX-SI_TT_E5_23W5C2LT4_S2_L005_R2_001.fastq.gz"
run_one "$ROOT/configs/e20_e24/e22b_gex.json" E22B \
  "$RAW/E22B_GEX/E22B_GEX-SI_TT_E6_23W5C2LT4_S3_L005_R1_001.fastq.gz" \
  "$RAW/E22B_GEX/E22B_GEX-SI_TT_E6_23W5C2LT4_S3_L005_R2_001.fastq.gz"
run_one "$ROOT/configs/e20_e24/e23tb_gex.json" E23TB \
  "$RAW/E23TB_GEX/E23TB_GEX-SI_TT_E7_23W5C2LT4_S4_L005_R1_001.fastq.gz" \
  "$RAW/E23TB_GEX/E23TB_GEX-SI_TT_E7_23W5C2LT4_S4_L005_R2_001.fastq.gz"
run_one "$ROOT/configs/e20_e24/e24t_gex.json" E24T \
  "$RAW/E24T_GEX/E24T_GEX-SI_TT_E8_23W5C2LT4_S5_L005_R1_001.fastq.gz" \
  "$RAW/E24T_GEX/E24T_GEX-SI_TT_E8_23W5C2LT4_S5_L005_R2_001.fastq.gz"
echo "ALL DONE $(date -Is)"
