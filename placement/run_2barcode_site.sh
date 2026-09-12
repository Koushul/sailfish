#!/usr/bin/env bash
# Build a here.now-ready 2-barcode localization site.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="${PYTHON:-python3}"

usage() {
  cat <<'EOF'
Usage:
  run_2barcode_site.sh --h5ad PATH --out DIR [options]

Required:
  --h5ad PATH              AnnData with obsm['ADT']
  --out DIR                Output site directory

Optional:
  --layout PATH            2D layout CSV (default: placement/layouts/layout_2d.csv)
  --feature-ref PATH       ADT feature CSV (default: refs/new_feature_ref_quant.csv)
  --dataset NAME           Dataset id in the UI (default: cells)
  --min-layout-umi N       Drop cells below this layout UMI (default: 10)
  --publish                Publish DIR with here.now after building
  --title TEXT
  --note TEXT

Example:
  placement/run_2barcode_site.sh \
      --h5ad /path/to/gex_adt.h5ad \
      --dataset MyExp \
      --out /tmp/myexp_site \
      --min-layout-umi 10
EOF
}

H5AD=""
OUT=""
EXTRA=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --h5ad) H5AD="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    *) EXTRA+=("$1"); shift ;;
  esac
done

if [[ -z "$H5AD" || -z "$OUT" ]]; then
  usage >&2
  exit 1
fi

exec "$PY" "$ROOT/placement/build_2barcode_site.py" --h5ad "$H5AD" --out "$OUT" "${EXTRA[@]}"
