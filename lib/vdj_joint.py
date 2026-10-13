"""One cell-level AnnData for a single animal: GEX of every library plus TCR/BCR.

All libraries come from the same mouse, so clone IDs are assigned across
libraries. BCR clonal families group heavy chains that share V, J and CDR3
length and are within 15% CDR3 nucleotide distance (single linkage), which
keeps somatically mutated members of one family together.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from .vdj import hamming
from .vdj_ocm import _as_int, tissue_population

FAMILY_MAX_DIST = 0.15
LIGHT = ("IGK", "IGL")


def _read_contigs(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        r["umis"] = _as_int(r.get("umis"))
    return rows


def _top(rows: list[dict], chains: tuple[str, ...]) -> list[dict]:
    return sorted((r for r in rows if r["chain"] in chains), key=lambda r: -r["umis"])


def _chain_fields(prefix: str, row: dict | None) -> dict:
    keys = ("cdr3", "cdr3_nt", "v_gene", "j_gene", "c_gene", "umis")
    if row is None:
        return {f"{prefix}_{k}": ("" if k != "umis" else 0) for k in keys}
    return {f"{prefix}_{k}": row.get(k, "") for k in keys}


def cell_table(library: str, receptor: str, contigs: list[dict]) -> dict[str, dict]:
    by_bc: dict[str, list[dict]] = defaultdict(list)
    for r in contigs:
        by_bc[r["barcode"]].append(r)
    cells = {}
    for bc, rows in by_bc.items():
        if receptor == "TCR":
            tra = _top(rows, ("TRA",))
            trb = _top(rows, ("TRB",))
            rec = {
                **_chain_fields("tra", tra[0] if tra else None),
                **_chain_fields("tra2", tra[1] if len(tra) > 1 else None),
                **_chain_fields("trb", trb[0] if trb else None),
                "tcr_paired": bool(tra and trb),
                "tcr_sample_id": rows[0].get("sample_id", ""),
            }
        else:
            igh = _top(rows, ("IGH",))
            light = _top(rows, LIGHT)
            rec = {
                **_chain_fields("igh", igh[0] if igh else None),
                **_chain_fields("light", light[0] if light else None),
                **_chain_fields("light2", light[1] if len(light) > 1 else None),
                "light_chain": light[0]["chain"] if light else "",
                "bcr_paired": bool(igh and light),
                "bcr_sample_id": rows[0].get("sample_id", ""),
            }
        cells[bc] = rec
    return cells


def _families(keys: list[tuple[str, str, str]]) -> dict[tuple[str, str, str], int]:
    """Single-linkage families over (V, J, CDR3 nt) heavy chains."""
    by_bin: dict[tuple[str, str, int], list[str]] = defaultdict(list)
    for v, j, nt in set(keys):
        by_bin[(v, j, len(nt))].append(nt)
    out: dict[tuple[str, str, str], int] = {}
    next_id = 1
    for (v, j, n), seqs in sorted(by_bin.items()):
        seqs = sorted(seqs)
        parent = list(range(len(seqs)))

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        limit = int(FAMILY_MAX_DIST * n)
        for a in range(len(seqs)):
            for b in range(a + 1, len(seqs)):
                if hamming(seqs[a], seqs[b]) <= limit:
                    parent[find(a)] = find(b)
        roots: dict[int, int] = {}
        for i, nt in enumerate(seqs):
            root = find(i)
            if root not in roots:
                roots[root] = next_id
                next_id += 1
            out[(v, j, nt)] = roots[root]
    return out


def _fill_vdj_columns(obs) -> None:
    for prefix in ("tra", "tra2", "trb", "igh", "light", "light2"):
        for k in ("cdr3", "cdr3_nt", "v_gene", "j_gene", "c_gene"):
            col = f"{prefix}_{k}"
            obs[col] = obs[col].fillna("").astype(str) if col in obs else ""
        ucol = f"{prefix}_umis"
        obs[ucol] = obs[ucol].fillna(0).astype(int) if ucol in obs else 0
    for col in ("tcr_paired", "bcr_paired"):
        obs[col] = obs[col].fillna(False).astype(bool) if col in obs else False
    for col in ("light_chain", "tcr_sample_id", "bcr_sample_id"):
        obs[col] = obs[col].fillna("").astype(str) if col in obs else ""
    for col in ("tcr_multiplet", "bcr_multiplet"):
        obs[col] = obs[col].fillna(False).astype(bool) if col in obs else False


def write_joint_h5ad(cfg: dict) -> Path:
    import anndata as ad
    import pandas as pd

    parts = []
    for lib in cfg["libraries"]:
        if not lib.get("gex_h5ad"):
            continue
        lib_id = lib["id"]
        gex = ad.read_h5ad(lib["gex_h5ad"])
        barcodes = [str(n).split("-")[0] for n in gex.obs_names]
        obs = gex.obs.copy()
        obs["library"] = lib_id
        obs["barcode"] = barcodes
        tissues = obs["sample_id"].astype(str).map(lambda s: tissue_population(s)[0])
        pops = obs["sample_id"].astype(str).map(lambda s: tissue_population(s)[1])
        obs["tissue"] = tissues.values
        obs["population"] = pops.values
        vdj_dir = Path(lib["vdj_dir"])
        for receptor in lib.get("receptors", []):
            path = vdj_dir / receptor / "filtered_contig_annotations.ocm.csv"
            if not path.exists():
                continue
            table = cell_table(lib_id, receptor.upper(), _read_contigs(path))
            frame = pd.DataFrame.from_dict(table, orient="index")
            frame = frame.reindex(barcodes)
            frame.index = obs.index
            for col in frame.columns:
                obs[col] = frame[col].values
            multiplets = vdj_dir / receptor / "vdj_multiplets.csv"
            if multiplets.exists():
                with multiplets.open(newline="") as fh:
                    flagged = {r["barcode"] for r in csv.DictReader(fh)}
                obs[f"{receptor.lower()}_multiplet"] = [bc in flagged for bc in barcodes]
        _fill_vdj_columns(obs)
        gex.obs = obs
        gex.obs_names = [f"{lib_id}|{bc}" for bc in barcodes]
        parts.append(gex)
    if not parts:
        raise ValueError("no libraries with gex_h5ad")
    joint = ad.concat(parts, join="inner", merge="same")
    obs = joint.obs
    obs["has_tcr"] = (obs["tra_cdr3_nt"] != "") | (obs["trb_cdr3_nt"] != "")
    obs["has_bcr"] = (obs["igh_cdr3_nt"] != "") | (obs["light_cdr3_nt"] != "")
    obs["receptor_conflict"] = obs["has_tcr"] & obs["has_bcr"]
    obs["vdj_multiplet"] = obs["tcr_multiplet"] | obs["bcr_multiplet"]
    obs["site"] = obs["library"].astype(str) + ":" + obs["sample_id"].astype(str)

    trb_key = obs["trb_v_gene"] + "|" + obs["trb_j_gene"] + "|" + obs["trb_cdr3_nt"]
    tcr_key = (
        obs["tra_v_gene"] + "|" + obs["tra_cdr3_nt"] + ";" + obs["tra2_v_gene"] + "|" + obs["tra2_cdr3_nt"]
        + ";" + trb_key
    )
    obs["trb_clone_id"] = ""
    obs.loc[obs["trb_cdr3_nt"] != "", "trb_clone_id"] = (
        "TRB" + trb_key[obs["trb_cdr3_nt"] != ""].astype("category").cat.codes.add(1).astype(str)
    )
    obs["tcr_clone_id"] = ""
    obs.loc[obs["has_tcr"], "tcr_clone_id"] = (
        "TCR" + tcr_key[obs["has_tcr"]].astype("category").cat.codes.add(1).astype(str)
    )
    igh_key = obs["igh_v_gene"] + "|" + obs["igh_j_gene"] + "|" + obs["igh_cdr3_nt"]
    bcr_key = igh_key + ";" + obs["light_v_gene"] + "|" + obs["light_cdr3_nt"]
    obs["bcr_clone_id"] = ""
    obs.loc[obs["has_bcr"], "bcr_clone_id"] = (
        "BCR" + bcr_key[obs["has_bcr"]].astype("category").cat.codes.add(1).astype(str)
    )
    has_igh = obs["igh_cdr3_nt"] != ""
    fam = _families(list(zip(obs.loc[has_igh, "igh_v_gene"], obs.loc[has_igh, "igh_j_gene"], obs.loc[has_igh, "igh_cdr3_nt"])))
    obs["bcr_family_id"] = ""
    obs.loc[has_igh, "bcr_family_id"] = [
        f"FAM{fam[(v, j, nt)]}"
        for v, j, nt in zip(obs.loc[has_igh, "igh_v_gene"], obs.loc[has_igh, "igh_j_gene"], obs.loc[has_igh, "igh_cdr3_nt"])
    ]
    for clone_col in ("trb_clone_id", "tcr_clone_id", "bcr_clone_id", "bcr_family_id"):
        mask = obs[clone_col] != ""
        sizes = obs.loc[mask].groupby(clone_col)[clone_col].transform("size")
        sites = obs.loc[mask].groupby(clone_col)["site"].transform("nunique")
        obs[f"{clone_col}_size"] = 0
        obs[f"{clone_col}_n_sites"] = 0
        obs.loc[mask, f"{clone_col}_size"] = sizes.astype(int)
        obs.loc[mask, f"{clone_col}_n_sites"] = sites.astype(int)
    for col in ("library", "sample_id", "tissue", "population", "site", "light_chain", "igh_c_gene", "trb_v_gene", "tra_v_gene", "igh_v_gene"):
        obs[col] = obs[col].astype(str).astype("category")
    joint.obs = obs
    joint.uns["organism"] = "Mus musculus"
    joint.uns["animal"] = cfg.get("animal", "")
    joint.uns["chemistry"] = "10xv3-5p"
    joint.uns["clone_rules"] = {
        "trb_clone_id": "TRB V|J|CDR3 nt, across all libraries of the animal",
        "tcr_clone_id": "exact TRA(1,2)+TRB V and CDR3 nt",
        "bcr_clone_id": "exact IGH V|J|CDR3 nt + light V|CDR3 nt",
        "bcr_family_id": f"IGH V, J, CDR3 length, single linkage at <= {FAMILY_MAX_DIST:.0%} CDR3 nt mismatch",
    }
    path = Path(cfg["joint_h5ad"])
    path.parent.mkdir(parents=True, exist_ok=True)
    joint.write_h5ad(path, compression="gzip")
    print(
        f"wrote {path} ({joint.n_obs} cells, {int(obs['has_tcr'].sum())} TCR, {int(obs['has_bcr'].sum())} BCR, "
        f"{int(obs['receptor_conflict'].sum())} with both)",
        flush=True,
    )
    return path


def cross_library_qc(cfg: dict) -> Path:
    """Match VDJ barcodes rejected by the GEX anchor against every library's GEX cells."""
    import anndata as ad

    gex = {}
    for lib in cfg["libraries"]:
        if lib.get("gex_h5ad"):
            names = ad.read_h5ad(lib["gex_h5ad"], backed="r").obs_names
            gex[lib["id"]] = {str(n).split("-")[0] for n in names}
    rows = []
    for lib in cfg["libraries"]:
        for receptor in lib.get("receptors", []):
            libdir = Path(lib["vdj_dir"]) / receptor
            for kind, fname in (("called", "filtered_contig_annotations.csv"), ("rejected", "vdj_cells_without_gex.csv")):
                path = libdir / fname
                if not path.exists():
                    continue
                with path.open(newline="") as fh:
                    barcodes = {r["barcode"] for r in csv.DictReader(fh)}
                by_oh: dict[str, set[str]] = defaultdict(set)
                for bc in barcodes:
                    by_oh[bc[7:9]].add(bc)
                for oh, bcs in sorted(by_oh.items()):
                    row = {"library": lib["id"], "receptor": receptor.upper(), "set": kind, "overhang": oh, "n": len(bcs)}
                    for other, cells in gex.items():
                        row[f"in_{other}_gex"] = len(bcs & cells)
                    rows.append(row)
    out = Path(cfg["qc_cross_library"])
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["library", "receptor", "set", "overhang", "n", *[f"in_{k}_gex" for k in gex]]
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}", flush=True)
    return out
