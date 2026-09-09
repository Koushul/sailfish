#!/usr/bin/env python3
"""48×48 microwell cell placement + entropy (parity with lucid-crystal site).

Algorithm used by https://lucid-crystal-kmqy.here.now/ :

  For each cell, each axis (row / column) is localized independently from two
  plates of 48 spatial-hash barcodes, then the MAP well is the product of the
  two axis choices.

  Per plate, per axis, barcode counts c ∈ ℕ^{48} with N = Σ c become a soft
  one-hot multinomial log-likelihood over the hypothesis "true index = k":

      p_hit  = e^β / (47 + e^β)
      p_miss = 1   / (47 + e^β)
      ll_k   = β · c_k − N · log(47 + e^β)     (default β = 3)

  Plates are combined with averaged empirical log-priors:

      lp_k = ll_k^(p1) + ll_k^(p2) + ½ (prior_k^(p1) + prior_k^(p2))

  Softmax(lp) → posterior π. Then:

      map            = argmax π                        (1..48)
      axis entropy   = −Σ π_k log2(π_k)               (bits; 0 if π_k≤1e-300)
      total entropy  = H(row) + H(col)                 (max 2·log2(48) ≈ 11.17)
      confidence     = max(π_row) · max(π_col)

  Spatial entropy (site "Weight by microwell proximity"): nearby alternate
  wells score lower than distant splits. The joint Π = π_row ⊗ π_col is
  not formed explicitly; the browser uses the separable identity
  Π ∗ K = (π_row ∗ g) ⊗ (π_col ∗ g) with a 1-D Gaussian g (σ = 1.5 wells):

      H_raw = −Σ π_row log2(π_row ∗ g) − Σ π_col log2(π_col ∗ g) + 2 log2(g(0))

  range-normalized so a uniform posterior still scores 2·log2(48).

  The live site stores per-plate ll / prior in cells.js and re-derives MAP,
  discrete entropy, and confidence in the browser with the same formulas
  (plate toggles drop one plate and use that plate's prior alone).
Examples
--------
  # Verify MAP/entropy parity against the published cells.js
  python cell_placement.py verify --cells-js https://lucid-crystal-kmqy.here.now/datasets/E28S/cells.js

  # Rebuild data.js from the 384-well oligo spreadsheet and check against the site
  python cell_placement.py layout-from-xlsx --xlsx layouts/chip_layout.xlsx --out layouts/data.js
  python cell_placement.py verify-layout --xlsx layouts/chip_layout.xlsx --data-js layouts/data.js

  # Assign from an AnnData with obsm['ADT'] + layout data.js
  python cell_placement.py from-h5ad --h5ad X.h5ad --data-js layouts/data.js --out assignments.csv
"""
from __future__ import annotations

import argparse
import base64
import json
import math
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np

GRID = 48
BETA = 3.0
MAX_ENTROPY_BITS = float(2.0 * math.log2(GRID))
DEFAULT_SPATIAL_SIGMA = 1.5  # wells; neighbor-scale Gaussian bandwidth
SEQ_SPLIT = "AGAATTCCA"
EPS = 1e-300
HERE = Path(__file__).resolve().parent
DEFAULT_CELLS_URL = "https://lucid-crystal-kmqy.here.now/datasets/E28S/cells.js"
DEFAULT_DATA_URL = "https://lucid-crystal-kmqy.here.now/data.js"
DEFAULT_LAYOUT_JS = HERE / "layouts" / "data.js"
DEFAULT_LAYOUT_XLSX = HERE / "layouts" / "chip_layout.xlsx"
DEFAULT_FEATURE_REF = HERE.parent / "refs" / "new_feature_ref_quant.csv"


# --------------------------------------------------------------------------- math
def axis_loglik(counts: np.ndarray, beta: float = BETA) -> np.ndarray:
    """Soft one-hot multinomial log-likelihood, shape (n, 48)."""
    counts = np.asarray(counts, dtype=np.float64)
    N = counts.sum(axis=1, keepdims=True)
    logZ = math.log(47.0 + math.exp(beta))
    return (beta * counts - N * logZ).astype(np.float32)


def axis_prior(counts: np.ndarray, alpha: float = 1.0) -> np.ndarray:
    """Empirical barcode log-prior from column sums, shape (48,)."""
    freq = np.asarray(counts, dtype=np.float64).sum(axis=0) + alpha
    freq /= freq.sum()
    return np.log(freq)


def combine_axis_log_post(
    ll_p1: np.ndarray,
    ll_p2: np.ndarray,
    prior_p1: np.ndarray,
    prior_p2: np.ndarray,
    use_p1: bool = True,
    use_p2: bool = True,
) -> np.ndarray:
    """Match site JS combineAxisLogPost (vectorized over cells)."""
    if not use_p1 and not use_p2:
        raise ValueError("Need at least one plate")
    ll_p1 = np.asarray(ll_p1, dtype=np.float64)
    ll_p2 = np.asarray(ll_p2, dtype=np.float64)
    prior_p1 = np.asarray(prior_p1, dtype=np.float64)
    prior_p2 = np.asarray(prior_p2, dtype=np.float64)
    lp = np.zeros_like(ll_p1, dtype=np.float64)
    prior = np.zeros(GRID, dtype=np.float64)
    n_plates = 0
    if use_p1:
        lp += ll_p1
        prior += prior_p1
        n_plates += 1
    if use_p2:
        lp += ll_p2
        prior += prior_p2
        n_plates += 1
    lp += prior / n_plates
    return lp


def normalize_log_post(lp: np.ndarray) -> np.ndarray:
    """Numerically stable softmax over the last axis (site normalizeLogPost)."""
    lp = np.asarray(lp, dtype=np.float64)
    m = lp.max(axis=-1, keepdims=True)
    p = np.exp(lp - m)
    p /= p.sum(axis=-1, keepdims=True)
    return p


def entropy_bits(p: np.ndarray) -> np.ndarray:
    """Shannon entropy in bits; mass ≤ 1e-300 contributes 0 (site entropyBits)."""
    p = np.asarray(p, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        logp = np.log(p) / math.log(2.0)
        logp[p <= EPS] = 0.0
        return -np.sum(p * logp, axis=-1)


def gaussian_kernel_1d(sigma: float) -> tuple[np.ndarray, int, float]:
    """Site gaussianKernel1D: integer-axis Gaussian, sum 1, k0 at center."""
    if sigma <= 0:
        raise ValueError("sigma must be > 0")
    rad = max(1, int(math.ceil(3.0 * sigma)))
    ax = np.arange(-rad, rad + 1, dtype=np.float64)
    k = np.exp(-(ax * ax) / (2.0 * sigma * sigma))
    k /= k.sum()
    return k, rad, float(k[rad])


def conv1d_same(p: np.ndarray, kern: np.ndarray, rad: int) -> np.ndarray:
    """Site conv1dSame: out[i] = Σ_j p[i-j] k[j+rad] for in-range i-j."""
    p = np.asarray(p, dtype=np.float64)
    n, g = p.shape
    out = np.zeros((n, g), dtype=np.float64)
    for j in range(-rad, rad + 1):
        w = kern[j + rad]
        if j >= 0:
            out[:, j:] += p[:, : g - j] * w
        else:
            out[:, : g + j] += p[:, -j:] * w
    return out


def spatial_entropy_bits_raw(
    row_post: np.ndarray,
    col_post: np.ndarray,
    sigma: float = DEFAULT_SPATIAL_SIGMA,
) -> np.ndarray:
    """Site spatialEntropyBitsRaw (separable 1-D Gaussian, not 2-D FFT)."""
    rp = np.asarray(row_post, dtype=np.float64)
    cp = np.asarray(col_post, dtype=np.float64)
    kern, rad, k0 = gaussian_kernel_1d(sigma)
    sr = conv1d_same(rp, kern, rad)
    sc = conv1d_same(cp, kern, rad)
    sr_sum = sr.sum(axis=1, keepdims=True)
    sc_sum = sc.sum(axis=1, keepdims=True)
    sr = np.where(sr_sum > EPS, sr / sr_sum, sr)
    sc = np.where(sc_sum > EPS, sc / sc_sum, sc)
    log2 = math.log(2.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        mask_r = (rp > EPS) & (sr > EPS)
        mask_c = (cp > EPS) & (sc > EPS)
        h = -np.sum(np.where(mask_r, rp * np.log(sr) / log2, 0.0), axis=1)
        h -= np.sum(np.where(mask_c, cp * np.log(sc) / log2, 0.0), axis=1)
    calib = 2.0 * math.log(k0) / log2
    return np.maximum(h + calib, 0.0)


def spatial_entropy_scale(sigma: float = DEFAULT_SPATIAL_SIGMA) -> float:
    """Multiplier mapping H_raw(uniform) → 2·log2(GRID) for this σ."""
    uni = np.full((1, GRID), 1.0 / GRID, dtype=np.float64)
    href = float(spatial_entropy_bits_raw(uni, uni, sigma=sigma)[0])
    if href <= EPS:
        return 1.0
    return MAX_ENTROPY_BITS / href


def spatial_entropy_bits(
    row_post: np.ndarray,
    col_post: np.ndarray,
    sigma: float = DEFAULT_SPATIAL_SIGMA,
) -> np.ndarray:
    """Range-normalized spatial entropy in bits (same ceiling as discrete H)."""
    raw = spatial_entropy_bits_raw(row_post, col_post, sigma=sigma)
    scaled = raw * spatial_entropy_scale(sigma)
    return np.minimum(scaled, MAX_ENTROPY_BITS)


def peak_norm_u8(p: np.ndarray) -> np.ndarray:
    peak = p.max(axis=-1, keepdims=True)
    inv = np.where(peak > EPS, 255.0 / peak, 0.0)
    return np.clip(np.round(p * inv), 0, 255).astype(np.uint8)


def place_from_loglik(
    row_ll_p1: np.ndarray,
    row_ll_p2: np.ndarray,
    col_ll_p1: np.ndarray,
    col_ll_p2: np.ndarray,
    row_prior_p1: np.ndarray,
    row_prior_p2: np.ndarray,
    col_prior_p1: np.ndarray,
    col_prior_p2: np.ndarray,
    use_p1: bool = True,
    use_p2: bool = True,
    spatial_sigma: float | None = DEFAULT_SPATIAL_SIGMA,
) -> dict[str, np.ndarray]:
    """Core placement: LLs + priors → MAP / entropy / confidence / posteriors."""
    rp = normalize_log_post(
        combine_axis_log_post(row_ll_p1, row_ll_p2, row_prior_p1, row_prior_p2, use_p1, use_p2)
    )
    cp = normalize_log_post(
        combine_axis_log_post(col_ll_p1, col_ll_p2, col_prior_p1, col_prior_p2, use_p1, use_p2)
    )
    row_h = entropy_bits(rp)
    col_h = entropy_bits(cp)
    out = {
        "map_row": rp.argmax(axis=1) + 1,
        "map_col": cp.argmax(axis=1) + 1,
        "confidence": rp.max(axis=1) * cp.max(axis=1),
        "row_entropy": row_h,
        "col_entropy": col_h,
        "total_entropy": row_h + col_h,
        "row_post": rp,
        "col_post": cp,
    }
    if spatial_sigma is not None:
        out["spatial_entropy"] = spatial_entropy_bits(rp, cp, sigma=spatial_sigma)
        out["spatial_sigma"] = np.asarray(spatial_sigma, dtype=np.float64)
    return out


def place_from_counts(
    row_p1: np.ndarray,
    row_p2: np.ndarray,
    col_p1: np.ndarray,
    col_p2: np.ndarray,
    beta: float = BETA,
    use_p1: bool = True,
    use_p2: bool = True,
    spatial_sigma: float | None = DEFAULT_SPATIAL_SIGMA,
) -> dict[str, np.ndarray]:
    """Counts → soft-one-hot LLs + empirical priors → placement."""
    row_ll_p1 = axis_loglik(row_p1, beta)
    row_ll_p2 = axis_loglik(row_p2, beta)
    col_ll_p1 = axis_loglik(col_p1, beta)
    col_ll_p2 = axis_loglik(col_p2, beta)
    out = place_from_loglik(
        row_ll_p1,
        row_ll_p2,
        col_ll_p1,
        col_ll_p2,
        axis_prior(row_p1),
        axis_prior(row_p2),
        axis_prior(col_p1),
        axis_prior(col_p2),
        use_p1=use_p1,
        use_p2=use_p2,
        spatial_sigma=spatial_sigma,
    )
    out["row_ll_p1"] = row_ll_p1
    out["row_ll_p2"] = row_ll_p2
    out["col_ll_p1"] = col_ll_p1
    out["col_ll_p2"] = col_ll_p2
    out["layout_umi"] = (
        np.asarray(row_p1) + np.asarray(row_p2) + np.asarray(col_p1) + np.asarray(col_p2)
    ).sum(axis=1).astype(np.int64)
    return out


# --------------------------------------------------------------------------- I/O
def _parse_js_payload(text: str, prefix: str) -> dict:
    text = text.strip()
    if not text.startswith(prefix):
        # allow bare JSON
        return json.loads(text.rstrip(";"))
    return json.loads(text.split("=", 1)[1].strip().rstrip(";"))


def _fetch_text(path_or_url: str) -> str:
    if path_or_url.startswith("http://") or path_or_url.startswith("https://"):
        req = urllib.request.Request(
            path_or_url,
            headers={"User-Agent": "cell_placement.py/1.0"},
        )
        with urllib.request.urlopen(req) as r:
            return r.read().decode("utf-8")
    return Path(path_or_url).read_text()


def load_cells_js(path_or_url: str) -> dict:
    return _parse_js_payload(_fetch_text(path_or_url), "window.CELL_DATA")


def load_data_js(path_or_url: str) -> dict:
    payload = _parse_js_payload(_fetch_text(path_or_url), "window.LAYOUT_DATA")
    by = {(b["plate"], b["axis"], int(b["index"])): b for b in payload["barcodes"]}
    return {"grid_size": int(payload["grid_size"]), "by": by, "barcodes": payload["barcodes"]}


def f32_from_b64(b64: str, n: int) -> np.ndarray:
    raw = base64.b64decode(b64)
    return np.frombuffer(raw, dtype=np.float32).reshape(n, GRID)


def u8_from_b64(b64: str, n: int) -> np.ndarray:
    raw = base64.b64decode(b64)
    return np.frombuffer(raw, dtype=np.uint8).reshape(n, GRID)


def extract_core(seq: str) -> str:
    if SEQ_SPLIT in seq:
        return seq.split(SEQ_SPLIT, 1)[1][:15]
    m = re.search(r"([ACGT]{15})B?A+$", seq)
    if not m:
        raise ValueError(f"Cannot parse oligo sequence: {seq[:40]}...")
    return m.group(1)


# --------------------------------------------------------------------------- verify
def verify_cells_js(cd: dict, rtol: float = 1e-4, atol: float = 1e-4) -> dict:
    """Recompute placement from embedded LLs/priors; compare to stored metrics."""
    n = int(cd["n_cells"])
    row_ll_p1 = f32_from_b64(cd["row_ll_p1_b64"], n)
    row_ll_p2 = f32_from_b64(cd["row_ll_p2_b64"], n)
    col_ll_p1 = f32_from_b64(cd["col_ll_p1_b64"], n)
    col_ll_p2 = f32_from_b64(cd["col_ll_p2_b64"], n)
    got = place_from_loglik(
        row_ll_p1,
        row_ll_p2,
        col_ll_p1,
        col_ll_p2,
        cd["row_prior_p1"],
        cd["row_prior_p2"],
        cd["col_prior_p1"],
        cd["col_prior_p2"],
    )

    ref_row = np.asarray(cd["map_row"], dtype=int)
    ref_col = np.asarray(cd["map_col"], dtype=int)
    ref_conf = np.asarray(cd["confidence"], dtype=np.float64)
    ref_rh = np.asarray(cd["row_entropy"], dtype=np.float64)
    ref_ch = np.asarray(cd["col_entropy"], dtype=np.float64)
    ref_th = np.asarray(cd["total_entropy"], dtype=np.float64)
    ref_rpost = u8_from_b64(cd["row_post_b64"], n)
    ref_cpost = u8_from_b64(cd["col_post_b64"], n)

    row_ok = bool(np.all(got["map_row"] == ref_row))
    col_ok = bool(np.all(got["map_col"] == ref_col))
    conf_err = float(np.max(np.abs(got["confidence"] - ref_conf)))
    rh_err = float(np.max(np.abs(got["row_entropy"] - ref_rh)))
    ch_err = float(np.max(np.abs(got["col_entropy"] - ref_ch)))
    th_err = float(np.max(np.abs(got["total_entropy"] - ref_th)))
    rpost_ok = bool(np.all(peak_norm_u8(got["row_post"]) == ref_rpost))
    cpost_ok = bool(np.all(peak_norm_u8(got["col_post"]) == ref_cpost))

    close = (
        np.allclose(got["confidence"], ref_conf, rtol=rtol, atol=atol)
        and np.allclose(got["row_entropy"], ref_rh, rtol=rtol, atol=atol)
        and np.allclose(got["col_entropy"], ref_ch, rtol=rtol, atol=atol)
        and np.allclose(got["total_entropy"], ref_th, rtol=rtol, atol=atol)
    )
    passed = row_ok and col_ok and rpost_ok and cpost_ok and close

    # single-plate modes (site plate toggles)
    plate_checks = {}
    for use_p1, use_p2, label in ((True, False, "p1_only"), (False, True, "p2_only")):
        g = place_from_loglik(
            row_ll_p1,
            row_ll_p2,
            col_ll_p1,
            col_ll_p2,
            cd["row_prior_p1"],
            cd["row_prior_p2"],
            cd["col_prior_p1"],
            cd["col_prior_p2"],
            use_p1=use_p1,
            use_p2=use_p2,
        )
        both = place_from_loglik(
            row_ll_p1,
            row_ll_p2,
            col_ll_p1,
            col_ll_p2,
            cd["row_prior_p1"],
            cd["row_prior_p2"],
            cd["col_prior_p1"],
            cd["col_prior_p2"],
        )
        plate_checks[label] = {
            "map_agree_with_both": float(
                np.mean((g["map_row"] == both["map_row"]) & (g["map_col"] == both["map_col"]))
            ),
            "median_total_entropy": float(np.median(g["total_entropy"])),
        }

    return {
        "passed": passed,
        "n_cells": n,
        "map_row_exact": row_ok,
        "map_col_exact": col_ok,
        "row_post_u8_exact": rpost_ok,
        "col_post_u8_exact": cpost_ok,
        "max_abs_confidence": conf_err,
        "max_abs_row_entropy": rh_err,
        "max_abs_col_entropy": ch_err,
        "max_abs_total_entropy": th_err,
        "max_entropy_bits_site": float(cd.get("max_entropy_bits", MAX_ENTROPY_BITS)),
        "max_entropy_bits_theory": MAX_ENTROPY_BITS,
        "plate_toggle_sanity": plate_checks,
    }


def _selftest() -> None:
    """Tiny synthetic check that counts → placement is deterministic."""
    rng = np.random.default_rng(0)
    n = 32
    row_p1 = np.zeros((n, GRID))
    row_p2 = np.zeros((n, GRID))
    col_p1 = np.zeros((n, GRID))
    col_p2 = np.zeros((n, GRID))
    true_r = rng.integers(0, GRID, size=n)
    true_c = rng.integers(0, GRID, size=n)
    for i in range(n):
        row_p1[i, true_r[i]] = 40
        row_p2[i, true_r[i]] = 35
        col_p1[i, true_c[i]] = 40
        col_p2[i, true_c[i]] = 35
        row_p1[i] += rng.integers(0, 2, GRID)
        row_p2[i] += rng.integers(0, 2, GRID)
        col_p1[i] += rng.integers(0, 2, GRID)
        col_p2[i] += rng.integers(0, 2, GRID)
    out = place_from_counts(row_p1, row_p2, col_p1, col_p2, beta=3.0)
    assert np.all(out["map_row"] == true_r + 1)
    assert np.all(out["map_col"] == true_c + 1)
    assert np.all(out["total_entropy"] < 0.5)
    assert np.all(out["confidence"] > 0.9)
    assert np.all(out["spatial_entropy"] < 0.5)
    # entropy of uniform is log2(48)
    uni = np.full(GRID, 1.0 / GRID)
    assert abs(float(entropy_bits(uni[None, :])[0]) - math.log2(GRID)) < 1e-12

    # Spatial entropy: same discrete H, nearby vs far ambiguity; range-normalized
    sigma = 1.5
    rp = np.zeros((2, GRID))
    cp = np.zeros((2, GRID))
    rp[0, 24] = 1.0
    cp[0, 10] = 0.5
    cp[0, 11] = 0.5
    rp[1, 24] = 1.0
    cp[1, 5] = 0.5
    cp[1, 40] = 0.5
    h_disc = entropy_bits(rp) + entropy_bits(cp)
    h_spat = spatial_entropy_bits(rp, cp, sigma=sigma)
    assert abs(h_disc[0] - h_disc[1]) < 1e-12
    assert h_spat[0] < h_spat[1] - 0.2, (h_spat[0], h_spat[1])
    uni = np.full((1, GRID), 1.0 / GRID)
    h_uni = float(spatial_entropy_bits(uni, uni, sigma=sigma)[0])
    assert abs(h_uni - MAX_ENTROPY_BITS) < 1e-6, h_uni
    rp0 = np.zeros((1, GRID)); rp0[0, 10] = 1.0
    cp0 = np.zeros((1, GRID)); cp0[0, 10] = 1.0
    assert abs(float(spatial_entropy_bits(rp0, cp0, sigma=sigma)[0])) < 1e-6
    print(
        f"selftest: ok  "
        f"(nearby H_sp={h_spat[0]:.3f} < far H_sp={h_spat[1]:.3f}; "
        f"discrete H={h_disc[0]:.3f}; uniform H_sp={h_uni:.3f})"
    )

# --------------------------------------------------------------------------- layout xlsx
_CHIP_LOC_RE = re.compile(r"^(ROW|COLUMN)\s+(\d+)$", re.I)
_SEQ_KEEP = set("ACGTacgtB")


def compact_oligo_sequence(seq: str) -> str:
    return "".join(ch for ch in str(seq) if ch in _SEQ_KEEP)


def layout_barcodes_from_xlsx(xlsx_path: str | Path) -> list[dict]:
    """Chip oligo spreadsheet → the 192 barcodes in site data.js form."""
    try:
        import pandas as pd
    except ImportError as e:
        raise SystemExit(f"layout-from-xlsx requires pandas/openpyxl: {e}") from e

    df = pd.read_excel(xlsx_path)
    df = df.dropna(subset=["Name", "Chip Location (Row/Column)"])
    barcodes = []
    for _, row in df.iterrows():
        loc = str(row["Chip Location (Row/Column)"]).strip()
        m = _CHIP_LOC_RE.match(loc)
        if not m:
            raise SystemExit(f"Unparseable chip location {loc!r} for {row['Name']}")
        axis = "row" if m.group(1).upper() == "ROW" else "column"
        index = int(m.group(2))
        plate_s = str(row["96 Source Plate #"]).strip()
        plate = int(plate_s.replace("Plate", "").strip())
        barcodes.append(
            {
                "name": str(row["Name"]),
                "plate": plate,
                "axis": axis,
                "index": index,
                "well": str(row["Well Position in 384 Source Plate"]),
                "sequence": compact_oligo_sequence(row["Sequence"]),
            }
        )
    if len(barcodes) != GRID * 4:
        raise SystemExit(f"Expected {GRID * 4} oligos, got {len(barcodes)}")
    return barcodes


def write_data_js(barcodes: list[dict], path: Path) -> None:
    payload = {"grid_size": GRID, "barcodes": barcodes}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("window.LAYOUT_DATA = " + json.dumps(payload, separators=(", ", ": ")) + ";\n")
    print(f"Wrote {path} ({len(barcodes)} barcodes)")


def verify_layout(xlsx_path: str, data_js: str) -> dict:
    layout = load_data_js(data_js)
    got = {(b["plate"], b["axis"], int(b["index"])): b for b in layout_barcodes_from_xlsx(xlsx_path)}
    ref = layout["by"]
    missing = [k for k in ref if k not in got]
    extra = [k for k in got if k not in ref]
    seq_mismatch = []
    name_mismatch = []
    well_mismatch = []
    for k, rb in ref.items():
        gb = got.get(k)
        if gb is None:
            continue
        if gb["sequence"] != rb["sequence"]:
            seq_mismatch.append(k)
        if gb["name"] != rb["name"]:
            name_mismatch.append(k)
        if gb["well"] != rb["well"]:
            well_mismatch.append(k)
    passed = not missing and not extra and not seq_mismatch and not name_mismatch and not well_mismatch
    return {
        "passed": passed,
        "n_xlsx": len(got),
        "n_data_js": len(ref),
        "missing_in_xlsx": missing,
        "extra_in_xlsx": extra,
        "seq_mismatch": seq_mismatch,
        "name_mismatch": name_mismatch,
        "well_mismatch": well_mismatch,
    }


# --------------------------------------------------------------------------- h5ad helper
def _adt_feature_names(adata, n_adt: int, feature_ref: str, adt_names: str | None) -> list[str]:
    import pandas as pd

    if adt_names:
        names = Path(adt_names).read_text().splitlines()
        if len(names) != n_adt:
            raise SystemExit(f"--adt-names has {len(names)} lines, ADT has {n_adt} columns")
        return names
    av = adata.uns.get("ADT_var")
    if av is not None:
        df = pd.DataFrame(av) if isinstance(av, dict) else av
        col = "feature_name" if "feature_name" in df.columns else df.columns[0]
        names = df[col].astype(str).tolist()
        if len(names) != n_adt:
            raise SystemExit(f"uns['ADT_var'] length {len(names)} != ADT columns {n_adt}")
        return names
    ref = pd.read_csv(feature_ref)
    if len(ref) != n_adt:
        raise SystemExit(
            f"h5ad has no uns['ADT_var'] and feature ref has {len(ref)} rows vs ADT {n_adt} columns; "
            "pass --adt-names (one feature name per ADT column)"
        )
    return ref["name"].astype(str).tolist()


def plate_axis_counts_from_h5ad(
    h5ad_path: str,
    data_js: str,
    feature_ref: str | None = None,
    adt_names: str | None = None,
):
    """Build (n × 48) count blocks ordered by data.js barcode indices."""
    try:
        import anndata as ad
        import pandas as pd
        from scipy import sparse
    except ImportError as e:
        raise SystemExit(f"from-h5ad requires anndata/pandas/scipy: {e}") from e

    if feature_ref is None:
        feature_ref = str(DEFAULT_FEATURE_REF)
    layout = load_data_js(data_js)
    adata = ad.read_h5ad(h5ad_path)
    if "ADT" not in adata.obsm:
        raise SystemExit("h5ad missing obsm['ADT']")
    ADT = adata.obsm["ADT"]
    ADT = ADT.toarray() if sparse.issparse(ADT) else np.asarray(ADT)
    names = _adt_feature_names(adata, ADT.shape[1], feature_ref, adt_names)
    ref = pd.read_csv(feature_ref)
    name_to_seq = dict(zip(ref["name"].astype(str), ref["sequence"].astype(str)))
    seq_to_col = {}
    for j, n in enumerate(names):
        seq = name_to_seq.get(n)
        if seq:
            seq_to_col[seq] = j

    blocks = {}
    for axis_js, axis_key in (("row", "row"), ("column", "col")):
        for plate in (1, 2):
            cols = []
            for idx in range(1, GRID + 1):
                bc = layout["by"][(plate, axis_js, idx)]
                core = extract_core(bc["sequence"])
                if core not in seq_to_col:
                    raise SystemExit(f"No ADT feature for {bc['name']} seq={core}")
                cols.append(seq_to_col[core])
            blocks[f"{axis_key}_p{plate}"] = ADT[:, cols].astype(np.float64)
    return blocks, adata.obs_names.astype(str).to_numpy()


def _write_csv(path: Path, names: list[str], out: dict) -> None:
    import csv

    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "barcode",
        "map_row",
        "map_col",
        "confidence",
        "row_entropy",
        "col_entropy",
        "total_entropy",
    ]
    if "spatial_entropy" in out:
        fields.append("spatial_entropy")
    if "layout_umi" in out:
        fields.append("layout_umi")
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for i, name in enumerate(names):
            row = {
                "barcode": name,
                "map_row": int(out["map_row"][i]),
                "map_col": int(out["map_col"][i]),
                "confidence": float(out["confidence"][i]),
                "row_entropy": float(out["row_entropy"][i]),
                "col_entropy": float(out["col_entropy"][i]),
                "total_entropy": float(out["total_entropy"][i]),
            }
            if "spatial_entropy" in out:
                row["spatial_entropy"] = float(out["spatial_entropy"][i])
            if "layout_umi" in out:
                row["layout_umi"] = int(out["layout_umi"][i])
            w.writerow(row)
    print(f"Wrote {path} (n={len(names)})")


# --------------------------------------------------------------------------- CLI
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_v = sub.add_parser("verify", help="Recompute from cells.js LLs and check parity")
    p_v.add_argument("--cells-js", default=DEFAULT_CELLS_URL)
    p_v.add_argument("--rtol", type=float, default=1e-4)
    p_v.add_argument("--atol", type=float, default=1e-4)

    p_lx = sub.add_parser("layout-from-xlsx", help="Write data.js from chip_layout.xlsx")
    p_lx.add_argument("--xlsx", default=str(DEFAULT_LAYOUT_XLSX))
    p_lx.add_argument("--out", default=str(DEFAULT_LAYOUT_JS))

    p_vl = sub.add_parser("verify-layout", help="Check xlsx oligos match data.js indices/sequences")
    p_vl.add_argument("--xlsx", default=str(DEFAULT_LAYOUT_XLSX))
    p_vl.add_argument("--data-js", default=str(DEFAULT_LAYOUT_JS))

    p_a = sub.add_parser("assign", help="Place cells from count matrices (.npy)")
    p_a.add_argument("--row-p1", required=True)
    p_a.add_argument("--row-p2", required=True)
    p_a.add_argument("--col-p1", required=True)
    p_a.add_argument("--col-p2", required=True)
    p_a.add_argument("--names", default=None, help="Text file of cell barcodes (one per line)")
    p_a.add_argument("--beta", type=float, default=BETA)
    p_a.add_argument("--spatial-sigma", type=float, default=DEFAULT_SPATIAL_SIGMA)
    p_a.add_argument(
        "--plates",
        choices=("both", "1", "2"),
        default="both",
        help="Which spatial-hash plates to use for localization (default: both)",
    )
    p_a.add_argument("--out", required=True)

    p_h = sub.add_parser("from-h5ad", help="Place cells from AnnData ADT + layout data.js")
    p_h.add_argument("--h5ad", required=True)
    p_h.add_argument("--data-js", default=str(DEFAULT_LAYOUT_JS))
    p_h.add_argument("--feature-ref", default=str(DEFAULT_FEATURE_REF))
    p_h.add_argument("--adt-names", default=None, help="Optional file of ADT column names (one per line)")
    p_h.add_argument("--beta", type=float, default=BETA)
    p_h.add_argument("--spatial-sigma", type=float, default=DEFAULT_SPATIAL_SIGMA)
    p_h.add_argument("--min-layout-umi", type=int, default=0)
    p_h.add_argument(
        "--plates",
        choices=("both", "1", "2"),
        default="both",
        help="Which spatial-hash plates to use for localization (default: both)",
    )
    p_h.add_argument("--out", required=True)

    p_s = sub.add_parser(
        "spatial-compare",
        help="Compare discrete vs spatial entropy on cells.js posteriors",
    )
    p_s.add_argument("--cells-js", default=DEFAULT_CELLS_URL)
    p_s.add_argument("--spatial-sigma", type=float, default=DEFAULT_SPATIAL_SIGMA)
    p_s.add_argument("--out", default=None, help="Optional CSV of per-cell entropies")

    sub.add_parser("selftest", help="Run synthetic deterministic checks")

    args = ap.parse_args(argv)

    if args.cmd == "selftest":
        _selftest()
        return 0

    if args.cmd == "layout-from-xlsx":
        write_data_js(layout_barcodes_from_xlsx(args.xlsx), Path(args.out))
        return 0

    if args.cmd == "verify-layout":
        report = verify_layout(args.xlsx, args.data_js)
        printable = {k: v for k, v in report.items() if k == "passed" or v}
        print(json.dumps(printable, indent=2, default=str))
        if not report["passed"]:
            print("LAYOUT MISMATCH", file=sys.stderr)
            return 1
        print("LAYOUT OK — xlsx oligos match data.js names, wells, sequences, and indices")
        return 0

    if args.cmd == "verify":
        print(f"Loading {args.cells_js}")
        cd = load_cells_js(args.cells_js)
        report = verify_cells_js(cd, rtol=args.rtol, atol=args.atol)
        print(json.dumps(report, indent=2))
        if not report["passed"]:
            print("PARITY FAILED", file=sys.stderr)
            return 1
        print("PARITY OK — MAP/posteriors exact; entropy/confidence within tolerance")
        return 0

    if args.cmd == "spatial-compare":
        print(f"Loading {args.cells_js}")
        cd = load_cells_js(args.cells_js)
        n = int(cd["n_cells"])
        got = place_from_loglik(
            f32_from_b64(cd["row_ll_p1_b64"], n),
            f32_from_b64(cd["row_ll_p2_b64"], n),
            f32_from_b64(cd["col_ll_p1_b64"], n),
            f32_from_b64(cd["col_ll_p2_b64"], n),
            cd["row_prior_p1"],
            cd["row_prior_p2"],
            cd["col_prior_p1"],
            cd["col_prior_p2"],
            spatial_sigma=args.spatial_sigma,
        )
        disc = got["total_entropy"]
        spat = got["spatial_entropy"]
        delta = disc - spat
        summary = {
            "n_cells": n,
            "spatial_sigma": args.spatial_sigma,
            "discrete_median": float(np.median(disc)),
            "spatial_median": float(np.median(spat)),
            "median_reduction_bits": float(np.median(delta)),
            "mean_reduction_bits": float(np.mean(delta)),
            "frac_spatial_lt_discrete": float(np.mean(spat < disc - 1e-9)),
            "frac_both_lt_1bit": {
                "discrete": float(np.mean(disc < 1.0)),
                "spatial": float(np.mean(spat < 1.0)),
            },
        }
        print(json.dumps(summary, indent=2))
        if args.out:
            _write_csv(Path(args.out), cd["obs_names"], got)
        return 0

    if args.cmd == "assign":
        row_p1 = np.load(args.row_p1)
        row_p2 = np.load(args.row_p2)
        col_p1 = np.load(args.col_p1)
        col_p2 = np.load(args.col_p2)
        n = row_p1.shape[0]
        if args.names:
            names = Path(args.names).read_text().splitlines()
            if len(names) != n:
                raise SystemExit(f"names length {len(names)} != n_cells {n}")
        else:
            names = [str(i) for i in range(n)]
        use_p1, use_p2 = {"both": (True, True), "1": (True, False), "2": (False, True)}[args.plates]
        out = place_from_counts(
            row_p1,
            row_p2,
            col_p1,
            col_p2,
            beta=args.beta,
            spatial_sigma=args.spatial_sigma,
            use_p1=use_p1,
            use_p2=use_p2,
        )
        if args.plates == "1":
            out["layout_umi"] = (np.asarray(row_p1) + np.asarray(col_p1)).sum(axis=1).astype(np.int64)
        elif args.plates == "2":
            out["layout_umi"] = (np.asarray(row_p2) + np.asarray(col_p2)).sum(axis=1).astype(np.int64)
        _write_csv(Path(args.out), names, out)
        te = out["total_entropy"]
        se = out["spatial_entropy"]
        print(
            f"discrete median={np.median(te):.3f}  spatial median={np.median(se):.3f}  "
            f"frac(H_sp<1)={(se < 1).mean():.3f}  frac(conf>0.9)={(out['confidence'] > 0.9).mean():.3f}"
        )
        return 0

    if args.cmd == "from-h5ad":
        blocks, names = plate_axis_counts_from_h5ad(
            args.h5ad, args.data_js, args.feature_ref, args.adt_names
        )
        if args.plates == "1":
            layout_umi = (blocks["row_p1"] + blocks["col_p1"]).sum(axis=1)
        elif args.plates == "2":
            layout_umi = (blocks["row_p2"] + blocks["col_p2"]).sum(axis=1)
        else:
            layout_umi = (
                blocks["row_p1"] + blocks["row_p2"] + blocks["col_p1"] + blocks["col_p2"]
            ).sum(axis=1)
        keep = layout_umi >= args.min_layout_umi
        blocks = {k: v[keep] for k, v in blocks.items()}
        names = names[keep].tolist()
        use_p1, use_p2 = {"both": (True, True), "1": (True, False), "2": (False, True)}[args.plates]
        out = place_from_counts(
            blocks["row_p1"],
            blocks["row_p2"],
            blocks["col_p1"],
            blocks["col_p2"],
            beta=args.beta,
            spatial_sigma=args.spatial_sigma,
            use_p1=use_p1,
            use_p2=use_p2,
        )
        if args.plates == "1":
            out["layout_umi"] = (blocks["row_p1"] + blocks["col_p1"]).sum(axis=1).astype(np.int64)
        elif args.plates == "2":
            out["layout_umi"] = (blocks["row_p2"] + blocks["col_p2"]).sum(axis=1).astype(np.int64)
        _write_csv(Path(args.out), names, out)
        te = out["total_entropy"]
        se = out["spatial_entropy"]
        print(
            f"plates={args.plates}  kept {len(names)} cells  discrete median={np.median(te):.3f}  "
            f"spatial median={np.median(se):.3f}  frac(conf>0.9)={(out['confidence'] > 0.9).mean():.3f}"
        )
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
