#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="${PYTHON:-python3}"
cd "$ROOT/placement"

echo "== selftest =="
"$PY" cell_placement.py selftest

echo "== verify-layout (xlsx vs data.js) =="
"$PY" cell_placement.py verify-layout \
  --xlsx layouts/chip_layout.xlsx \
  --data-js layouts/data.js

echo "== verify algorithm vs lucid-crystal E28S cells.js =="
"$PY" cell_placement.py verify \
  --cells-js https://lucid-crystal-kmqy.here.now/datasets/E28S/cells.js
