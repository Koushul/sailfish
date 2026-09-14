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


def _obs_vector(obs: h5py.Group, name: str) -> np.ndarray:
    g = obs[name]
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


def load_placed(path: str) -> PlacedCounts:
    targets = core_targets()
    mouse = [t.mouse for t in targets]
    with h5py.File(path, "r") as f:
        n_obs = int(f["obs"]["_index"].shape[0])
        genes = np.array([_decode(x) for x in f["var"]["_index"][:]])
        gix = {g: i for i, g in enumerate(genes)}
        missing = [g for g in mouse if g not in gix]
        if missing:
            raise KeyError(f"HIF panel genes missing from var: {missing}")
        panel_idx = np.array([gix[g] for g in mouse], dtype=np.int32)
        s_idx = np.array([gix[g] for g in S_GENES if g in gix], dtype=np.int32)
        g2m_idx = np.array([gix[g] for g in G2M_GENES if g in gix], dtype=np.int32)
        spliced = _csr_layer(f["layers"]["spliced"], n_obs, genes.size)
        unspliced = _csr_layer(f["layers"]["unspliced"], n_obs, genes.size)
        obs = f["obs"]
        cell_id = _obs_vector(obs, "_index")
        sample = _obs_vector(obs, "sample")
        cell_group = _obs_vector(obs, "cell_group")
        hist_state = (
            _obs_vector(obs, "hypoxia_kinetics_state")
            if "hypoxia_kinetics_state" in obs
            else np.array(["NA"] * n_obs, dtype=object)
        )
        hist_theta = (
            np.asarray(obs["theta_normoxic"][:], dtype=np.float64)
            if "theta_normoxic" in obs
            else np.full(n_obs, np.nan)
        )

    L = np.asarray(spliced.sum(axis=1)).ravel()
    S = spliced[:, panel_idx].toarray()
    U = unspliced[:, panel_idx].toarray()
    scale = np.median(L) / np.clip(L, 1.0, None)
    log_s = np.log1p(spliced[:, s_idx].toarray() * scale[:, None]) if s_idx.size else np.zeros((n_obs, 1))
    log_g = np.log1p(spliced[:, g2m_idx].toarray() * scale[:, None]) if g2m_idx.size else np.zeros((n_obs, 1))
    cycle_s = log_s.mean(axis=1)
    cycle_g2m = log_g.mean(axis=1)
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
    )
