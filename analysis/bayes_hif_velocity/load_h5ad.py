#!/usr/bin/env python3
"""Read spliced/unspliced counts from a placed AnnData file (read-only)."""
from __future__ import annotations

from dataclasses import dataclass

import h5py
import numpy as np
from scipy.sparse import csr_matrix

from hif_targets import core_targets
from tirosh_mouse import G2M_GENES, S_GENES


def _decode(x) -> str:
    return x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x)


def _ann_vector(group: h5py.Group, name: str) -> np.ndarray:
    g = group[name]
    if isinstance(g, h5py.Group) and "categories" in g and "codes" in g:
        cats = np.array([_decode(x) for x in g["categories"][:]])
        codes = np.asarray(g["codes"][:], dtype=np.int64)
        out = np.empty(codes.shape[0], dtype=object)
        ok = codes >= 0
        out[ok] = cats[codes[ok]]
        out[~ok] = "NA"
        return out
    a = g[:]
    if a.dtype.kind in ("S", "O", "U"):
        return np.array([_decode(x) for x in a], dtype=object)
    return np.asarray(a)


def _obs_vector(obs: h5py.Group, name: str) -> np.ndarray:
    return _ann_vector(obs, name)


def _csr_layer(group: h5py.Group, n_obs: int, n_var: int) -> csr_matrix:
    data = np.asarray(group["data"][:], dtype=np.float64)
    indices = np.asarray(group["indices"][:], dtype=np.int32)
    indptr = np.asarray(group["indptr"][:], dtype=np.int32)
    return csr_matrix((data, indices, indptr), shape=(n_obs, n_var))


@dataclass
class PlacedCounts:
    cell_id: np.ndarray
    sample: np.ndarray
    cell_group: np.ndarray
    L: np.ndarray
    genes: np.ndarray
    S: np.ndarray
    U: np.ndarray
    cycle_s: np.ndarray
    cycle_g2m: np.ndarray
    n_s_genes: int
    n_g2m_genes: int
    historical_state: np.ndarray
    historical_theta: np.ndarray
    missing_panel: tuple[str, ...]


def _symbol_index(var: h5py.Group, gene_symbol_col: str | None) -> dict[str, int]:
    n_var = int(var["_index"].shape[0])
    primary = np.array([_decode(x) for x in var["_index"][:]])
    gix = {g: i for i, g in enumerate(primary)}
    if gene_symbol_col and gene_symbol_col in var:
        alt = _ann_vector(var, gene_symbol_col)
        for i, g in enumerate(alt):
            if g and g != "NA" and g not in gix:
                gix[g] = i
        for i, g in enumerate(primary):
            gix.setdefault(g, i)
    elif n_var:
        pass
    return gix


def load_counts(
    path: str,
    *,
    sample_col: str = "sample",
    lineage_col: str = "cell_group",
    gene_symbol_col: str | None = None,
    historical_state_col: str | None = "hypoxia_kinetics_state",
    historical_theta_col: str | None = "theta_normoxic",
    require_all_panel: bool = False,
) -> PlacedCounts:
    targets = core_targets()
    mouse = [t.mouse for t in targets]
    with h5py.File(path, "r") as f:
        n_obs = int(f["obs"]["_index"].shape[0])
        var = f["var"]
        gix = _symbol_index(var, gene_symbol_col)
        n_var = int(var["_index"].shape[0])
        missing = [g for g in mouse if g not in gix]
        if missing and require_all_panel:
            raise KeyError(f"HIF panel genes missing from var: {missing}")
        panel_idx = np.array([gix[g] if g in gix else -1 for g in mouse], dtype=np.int32)
        s_idx = np.array([gix[g] for g in S_GENES if g in gix], dtype=np.int32)
        g2m_idx = np.array([gix[g] for g in G2M_GENES if g in gix], dtype=np.int32)
        spliced = _csr_layer(f["layers"]["spliced"], n_obs, n_var)
        unspliced = _csr_layer(f["layers"]["unspliced"], n_obs, n_var)
        obs = f["obs"]
        cell_id = _obs_vector(obs, "_index")
        if sample_col not in obs:
            raise KeyError(f"obs is missing sample column {sample_col!r}")
        sample = _obs_vector(obs, sample_col)
        if lineage_col and lineage_col in obs:
            cell_group = _obs_vector(obs, lineage_col)
        else:
            cell_group = np.array(["unlabeled"] * n_obs, dtype=object)
        hist_col = historical_state_col if historical_state_col and historical_state_col in obs else None
        hist_state = (
            _obs_vector(obs, hist_col) if hist_col else np.array(["NA"] * n_obs, dtype=object)
        )
        if historical_theta_col and historical_theta_col in obs:
            hist_theta = np.asarray(obs[historical_theta_col][:], dtype=np.float64)
        else:
            hist_theta = np.full(n_obs, np.nan)

    L = np.asarray(spliced.sum(axis=1)).ravel()
    S = np.zeros((n_obs, len(mouse)), dtype=np.float64)
    U = np.zeros((n_obs, len(mouse)), dtype=np.float64)
    present = panel_idx >= 0
    if np.any(present):
        S[:, present] = spliced[:, panel_idx[present]].toarray()
        U[:, present] = unspliced[:, panel_idx[present]].toarray()
    cpm_scale = 1e4 / np.clip(L, 1.0, None)
    log_s = log_g = None
    log_s = spliced[:, s_idx].toarray() * cpm_scale[:, None] if s_idx.size else np.zeros((n_obs, 1))
    log_g = spliced[:, g2m_idx].toarray() * cpm_scale[:, None] if g2m_idx.size else np.zeros((n_obs, 1))
    cycle_s = np.log1p(log_s).mean(axis=1)
    cycle_g2m = np.log1p(log_g).mean(axis=1)
    return PlacedCounts(
        cell_id=cell_id,
        sample=sample,
        cell_group=cell_group,
        L=L,
        genes=np.array(mouse, dtype=object),
        S=S,
        U=U,
        cycle_s=cycle_s,
        cycle_g2m=cycle_g2m,
        n_s_genes=int(s_idx.size),
        n_g2m_genes=int(g2m_idx.size),
        historical_state=hist_state,
        historical_theta=hist_theta,
        missing_panel=tuple(missing),
    )


def load_placed(path: str) -> PlacedCounts:
    return load_counts(
        path,
        sample_col="sample",
        lineage_col="cell_group",
        gene_symbol_col=None,
        historical_state_col="hypoxia_kinetics_state",
        historical_theta_col="theta_normoxic",
        require_all_panel=True,
    )
