#!/bin/bash
# Call TCR/BCR for E32–E36 and write the joint GEX + VDJ object.
# Requires the OCM-split GEX h5ad from run_e32_e36_gex.sh.

set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
PY="${PY:-/ix1/ylee/kor11/tools/af_tutorial/conda_env/bin/python}"
LOGS=/ix1/ylee/kor11/TCR_BCR/logs
mkdir -p "$LOGS" /ix1/ylee/kor11/TCR_BCR/runs
cd "$ROOT"
EXTRA=("$@")

run_one() {
  local cfg="$1" receptor="$2"
  local name
  name="$(basename "$cfg" .json)_${receptor}"
  echo "START $(date -Is) $name"
  "$PY" "$ROOT/run.py" "$cfg" --vdj-receptor "$receptor" "${EXTRA[@]}" > "$LOGS/${name}.log" 2>&1
  echo "DONE $(date -Is) $name"
}

C="$ROOT/configs/e32_e36"
run_one "$C/e32t.json" tcr &
run_one "$C/e33t.json" tcr &
run_one "$C/e34b.json" bcr &
run_one "$C/e35tb.json" tcr &
run_one "$C/e35tb.json" bcr &
run_one "$C/e36t.json" tcr &
wait
echo "VDJ DONE $(date -Is)"
"$PY" "$ROOT/run.py" "$C/ocm_split.json" --split-vdj-ocm > "$LOGS/e32_e36_ocm_split.log" 2>&1
echo "SPLIT DONE $(date -Is)"
