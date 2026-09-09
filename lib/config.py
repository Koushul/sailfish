from __future__ import annotations

from pathlib import Path

from .demux_ocm import GEMX_OCM_DEFAULT

GEMX_OVERHANG_TO_OB = {"GT": "OB1", "CA": "OB2", "TC": "OB3", "AG": "OB4"}
GEMX_OB_TO_OVERHANG = {v: k for k, v in GEMX_OVERHANG_TO_OB.items()}


def _p(v):
    return None if v is None else Path(v)


def parse_ocm_samples(ocm: dict | None) -> list[dict]:
    ocm = ocm or {}
    raw = ocm.get("samples")
    if not raw:
        return list(GEMX_OCM_DEFAULT)
    if isinstance(raw, dict):
        out = []
        for key, val in raw.items():
            key = str(key)
            if isinstance(val, dict):
                sample_id = val.get("sample_id", key)
                overhang = val.get("overhang", GEMX_OB_TO_OVERHANG.get(key, key))
                desc = val.get("description", sample_id)
                ob = val.get("ocm_barcode_id", GEMX_OVERHANG_TO_OB.get(overhang, key))
            else:
                sample_id = str(val)
                if key in GEMX_OVERHANG_TO_OB:
                    overhang, ob = key, GEMX_OVERHANG_TO_OB[key]
                elif key in GEMX_OB_TO_OVERHANG:
                    ob, overhang = key, GEMX_OB_TO_OVERHANG[key]
                else:
                    overhang, ob = key, key
                desc = sample_id
            out.append(
                {
                    "ocm_barcode_id": ob,
                    "sample_id": sample_id,
                    "description": desc,
                    "overhang": overhang,
                }
            )
        return out
    out = []
    for s in raw:
        overhang = s["overhang"]
        ob = s.get("ocm_barcode_id") or GEMX_OVERHANG_TO_OB.get(overhang, s.get("sample_id"))
        sample_id = s.get("sample_id", ob)
        out.append(
            {
                "ocm_barcode_id": ob,
                "sample_id": sample_id,
                "description": s.get("description", sample_id),
                "overhang": overhang,
            }
        )
    return out


def normalize_config(cfg: dict) -> dict:
    cfg = dict(cfg)
    if "output" not in cfg:
        raise ValueError("config needs 'output'")
    if "gex" not in cfg:
        raise ValueError("config needs 'gex' (reads1, reads2, index, chemistry)")
    cfg.setdefault("sample", Path(cfg["output"]).name)
    cfg.setdefault("mode", "quant")
    cfg.setdefault("threads", 16)
    cfg.setdefault("skip_quant", False)
    if "min_gex_umi" not in cfg:
        cfg["min_gex_umi"] = int((cfg.get("filter") or {}).get("min_gex_umi", 500))
    cfg.setdefault("alevin_fry_home", str(Path(cfg["output"]) / "af_home"))
    gex = dict(cfg["gex"])
    gex.setdefault("min_reads", 10)
    gex.setdefault("resolution", "cr-like")
    gex.setdefault("chemistry", "10xv3")
    cfg["gex"] = gex
    if cfg.get("adt"):
        adt = dict(cfg["adt"])
        adt.setdefault("min_reads", 10)
        adt.setdefault("resolution", "cr-like")
        cfg["adt"] = adt
    if cfg["mode"] == "ocm":
        ocm = dict(cfg.get("ocm") or {})
        ocm.setdefault("overhang_start", 7)
        ocm.setdefault("overhang_len", 2)
        ocm.setdefault("include_unassigned", False)
        ocm["samples"] = parse_ocm_samples(ocm)
        cfg["ocm"] = ocm
    return cfg
