#!/bin/bash
# Quantify E32–E36 GEM-X 5' GEX and split each library by OCM well.
# Second lupus mouse. Well order matches E20–E24.

set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
PY="${PY:-/ix1/ylee/kor11/tools/af_tutorial/conda_env/bin/python}"
STAGE="${STAGE:-/scratch/slurm-4343839/e32_e36_gex}"
INDEX_SRC="/scratch/slurm-4343839/e20_e24_gex/index"
if [[ ! -s "$INDEX_SRC/t2g_3col.tsv" ]]; then
  INDEX_SRC="/ix1/ylee/shared/E18S_E19S_E37S1_E37S2/E37S2_sailfish/refs/mouse-2024-A_splici/index"
fi
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
  if [[ -s "$DEST/$name/gex_quant/af_quant/alevin/quants.h5ad" && -s "$DEST/$name/ocm/${name}_gex_adt_ocm.h5ad" ]]; then
    echo "SKIP $(date -Is) $name"
    return 0
  fi
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
"$PY" - << PY
import gzip
from pathlib import Path
src_r1 = Path("$RAW/E32T_GEX/E32T_GEX-SI_TT_G1_23W5C2LT4_S6_L005_R1_001.fastq.gz")
src_r2 = Path("$RAW/E32T_GEX/E32T_GEX-SI_TT_G1_23W5C2LT4_S6_L005_R2_001.fastq.gz")
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
staged="$(write_staged_config "$ROOT/configs/e32_e36/e32t_gex.json" pilot "$pilot_out")"
"$PY" "$ROOT/run.py" "$staged" > "$LOGDIR/e32_pilot_gex.log" 2>&1
"$PY" - << PY
import json
from pathlib import Path
info = json.loads(Path("$pilot_out/gex_quant/af_map/map_info.json").read_text())
pct = float(info["percent_mapped"])
print(f"pilot mapped {info['num_mapped']}/{info['num_reads']} ({pct}%)", flush=True)
cmd = info.get("cmdline", "")
if "chromium_v3_5p" not in cmd:
    raise SystemExit(f"expected chromium_v3_5p in piscem command: {cmd}")
if pct < 50:
    raise SystemExit(f"mapping rate {pct} is below 50")
PY
rm -f "$STAGE/fastq/pilot_R1.fastq.gz" "$STAGE/fastq/pilot_R2.fastq.gz"
rm -rf "$pilot_out/gex_quant/af_map"
echo "PILOT OK $(date -Is)"

run_one "$ROOT/configs/e32_e36/e32t_gex.json" E32T \
  "$RAW/E32T_GEX/E32T_GEX-SI_TT_G1_23W5C2LT4_S6_L005_R1_001.fastq.gz" \
  "$RAW/E32T_GEX/E32T_GEX-SI_TT_G1_23W5C2LT4_S6_L005_R2_001.fastq.gz"
run_one "$ROOT/configs/e32_e36/e33t_gex.json" E33T \
  "$RAW/E33T_GEX/E33T_GEX-SI_TT_G2_23W5C2LT4_S7_L005_R1_001.fastq.gz" \
  "$RAW/E33T_GEX/E33T_GEX-SI_TT_G2_23W5C2LT4_S7_L005_R2_001.fastq.gz"
run_one "$ROOT/configs/e32_e36/e34b_gex.json" E34B \
  "$RAW/E34B_GEX/E34B_GEX-SI_TT_G3_23W5C2LT4_S1_L006_R1_001.fastq.gz" \
  "$RAW/E34B_GEX/E34B_GEX-SI_TT_G3_23W5C2LT4_S1_L006_R2_001.fastq.gz"
run_one "$ROOT/configs/e32_e36/e35tb_gex.json" E35TB \
  "$RAW/E35TB_GEX/E35TB_GEX-SI_TT_G4_23W5C2LT4_S2_L006_R1_001.fastq.gz" \
  "$RAW/E35TB_GEX/E35TB_GEX-SI_TT_G4_23W5C2LT4_S2_L006_R2_001.fastq.gz"
run_one "$ROOT/configs/e32_e36/e36t_gex.json" E36T \
  "$RAW/E36T_GEX/E36T_GEX-SI_TT_G5_23W5C2LT4_S3_L006_R1_001.fastq.gz" \
  "$RAW/E36T_GEX/E36T_GEX-SI_TT_G5_23W5C2LT4_S3_L006_R2_001.fastq.gz"
echo "GEX DONE $(date -Is)"
bash "$ROOT/run_e32_e36_vdj.sh"
echo "ALL DONE $(date -Is)"
