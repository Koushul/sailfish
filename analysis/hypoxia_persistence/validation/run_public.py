#!/usr/bin/env python3
"""Validate the Tumor HIF module on public labeled datasets (priority order)."""
from __future__ import annotations

import gzip
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.metrics import roc_auc_score

from score import (
    gene_auroc,
    human_theta_genes,
    log1p_cp10k,
    module_score,
    summarize_split,
)

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT = HERE / "results"
GPL6884 = DATA / "GPL6884.annot.gz"
SERIES30019 = DATA / "GSE30019_series_matrix.txt.gz"


def load_theta_panel():
    table = human_theta_genes()
    return table["human"].tolist(), table["weight"].to_numpy(dtype=np.float64), table


def extract_gene_matrix_tsv(path: Path, genes: list[str]) -> tuple[np.ndarray, np.ndarray]:
    want = {g.upper() for g in genes}
    found = {}
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cells = header[1:]
        for line in f:
            parts = line.rstrip("\n").split("\t", 1)
            if not parts or parts[0].upper() not in want:
                continue
            vals = np.fromstring(parts[1], sep="\t", dtype=np.float64) if len(parts) > 1 else np.array([])
            found[parts[0].upper()] = vals
    X = np.column_stack([found[g.upper()] if g.upper() in found else np.zeros(len(cells)) for g in genes])
    return X, np.array(cells)


def gse292771() -> pd.DataFrame:
    genes, weights, _ = load_theta_panel()
    files = {
        "WT_Nx": DATA / "GSE292771" / "GSM8864979_1-48Nx_20230427202845_TCM.tsv.gz",
        "WT_Hx": DATA / "GSE292771" / "GSM8864980_2-48Hx_20230427205149_TCM.tsv.gz",
        "HIF1KO_Nx": DATA / "GSE292771" / "GSM8864977_1-KO_HIF-1a_Nx_20221220103821_TCM.tsv.gz",
        "HIF1KO_Hx": DATA / "GSE292771" / "GSM8864978_2-KO_HIF-1a_Hx_20221220110433_TCM.tsv.gz",
        "HIF2KO_Nx": DATA / "GSE292771" / "GSM8864975_1-KO_HIF-2a_Nx_20221220111305_TCM.tsv.gz",
        "HIF2KO_Hx": DATA / "GSE292771" / "GSM8864976_2-KO_HIF-2a_Hx_20221220112633_TCM.tsv.gz",
    }
    blocks = []
    for lab, path in files.items():
        X, cells = extract_gene_matrix_tsv(path, genes)
        logx = log1p_cp10k(X)
        genotype, ox = lab.split("_")
        blocks.append(
            pd.DataFrame(
                {
                    "cell": [f"{lab}:{c}" for c in cells],
                    "genotype": genotype,
                    "oxygen": ox,
                    "dataset": "GSE292771",
                    **{g: logx[:, i] for i, g in enumerate(genes)},
                }
            )
        )
    df = pd.concat(blocks, ignore_index=True)
    X = df[genes].to_numpy(dtype=np.float64)
    df["hif_module"] = module_score(X, weights)
    rows = []
    for geno in ["WT", "HIF1KO", "HIF2KO"]:
        sub = df[df.genotype == geno]
        met = summarize_split(sub.hif_module.to_numpy(), sub.oxygen.to_numpy(), "Hx", "Nx")
        met["dataset"] = "GSE292771"
        met["split"] = f"{geno} Hx vs Nx (module)"
        rows.append(met)
        y = (sub.oxygen == "Hx").to_numpy()
        gtab = gene_auroc(sub[genes].to_numpy(), y, genes)
        gtab["dataset"] = "GSE292771"
        gtab["split"] = f"{geno} Hx vs Nx"
        gtab.to_csv(OUT / f"GSE292771_gene_auroc_{geno}.csv", index=False)
    # module should beat median single-gene AUROC in WT
    wt = pd.read_csv(OUT / "GSE292771_gene_auroc_WT.csv")
    rows.append(
        {
            "dataset": "GSE292771",
            "split": "WT module vs median gene AUROC",
            "auroc": rows[0]["auroc"],
            "median_gene_auroc": float(wt.auroc.median()),
            "best_gene": wt.sort_values("auroc").iloc[-1]["gene"],
            "best_gene_auroc": float(wt.auroc.max()),
            "n": rows[0]["n"],
            "n_pos": rows[0]["n_pos"],
        }
    )
    return pd.DataFrame(rows), df


def ensure_gpl6884() -> pd.DataFrame:
    if not GPL6884.exists():
        import urllib.request

        url = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6884/annot/GPL6884.annot.gz"
        urllib.request.urlretrieve(url, GPL6884)
    rows = []
    with gzip.open(GPL6884, "rt", errors="replace") as f:
        in_table = False
        header = None
        for line in f:
            if line.startswith("ID\t") or line.startswith("ID,"):
                in_table = True
                header = line.rstrip("\n").split("\t")
                continue
            if not in_table or line.startswith("#") or line.startswith("^"):
                continue
            parts = line.rstrip("\n").split("\t")
            rec = dict(zip(header, parts))
            rows.append(rec)
            if len(rows) > 60000:
                break
    ann = pd.DataFrame(rows)
    symbol_col = next((c for c in ann.columns if c.lower() in {"gene symbol", "genesymbol", "symbol"}), None)
    if symbol_col is None:
        symbol_col = next((c for c in ann.columns if "symbol" in c.lower()), ann.columns[1])
    return ann.rename(columns={header[0] if False else ann.columns[0]: "ID", symbol_col: "symbol"})[["ID", "symbol"]]


def gse30019() -> pd.DataFrame:
    genes, weights, _ = load_theta_panel()
    if not SERIES30019.exists():
        import urllib.request

        urllib.request.urlretrieve(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30019/matrix/GSE30019_series_matrix.txt.gz",
            SERIES30019,
        )
    mat = pd.read_csv(SERIES30019, sep="\t", comment="!", index_col=0)
    mat = mat.dropna(how="all")
    samples = [
        ("GSM742971", 0),
        ("GSM742972", 0),
        ("GSM742973", 0),
        ("GSM742974", 1),
        ("GSM742975", 1),
        ("GSM742976", 1),
        ("GSM742977", 4),
        ("GSM742978", 4),
        ("GSM742979", 4),
        ("GSM742980", 8),
        ("GSM742981", 8),
        ("GSM742982", 8),
        ("GSM742983", 12),
        ("GSM742984", 12),
        ("GSM742985", 12),
        ("GSM742986", 24),
        ("GSM742987", 24),
        ("GSM742988", 24),
    ]
    cols = [s for s, _ in samples]
    mat = mat[cols].apply(pd.to_numeric, errors="coerce")
    ann = ensure_gpl6884()
    ann["symbol"] = ann["symbol"].astype(str).str.split("///").str[0].str.strip().str.upper()
    want = {g.upper() for g in genes}
    ann = ann[ann["symbol"].isin(want)]
    expr = mat.join(ann.set_index("ID")["symbol"], how="inner")
    expr = expr.groupby("symbol").mean()
    X = np.vstack([expr.reindex([g.upper()]).to_numpy() for g in genes]).T  # samples x genes
    X = np.nan_to_num(X, nan=np.nanmedian(X))
    # z within this bulk set
    score = module_score(np.log2(np.clip(X, 1, None)), weights)
    hours = np.array([h for _, h in samples], dtype=float)
    rho = float(stats.spearmanr(hours, score).correlation)
    # 0h (still hypoxic) vs 24h reox
    lab = np.where(hours == 0, "Hx", np.where(hours == 24, "Nx", "mid"))
    met = summarize_split(score, lab, "Hx", "Nx")
    met.update({"dataset": "GSE30019", "split": "0h hypoxia vs 24h reox (bulk module)", "spearman_hours": rho})
    per = pd.DataFrame({"sample": cols, "hours": hours, "hif_module": score})
    per.to_csv(OUT / "GSE30019_sample_scores.csv", index=False)
    gene_change = []
    for i, g in enumerate(genes):
        gene_change.append(
            {
                "gene": g,
                "mean_0h": float(X[hours == 0, i].mean()),
                "mean_24h": float(X[hours == 24, i].mean()),
                "log2fc_24_vs_0": float(np.log2(X[hours == 24, i].mean() + 1) - np.log2(X[hours == 0, i].mean() + 1)),
            }
        )
    pd.DataFrame(gene_change).to_csv(OUT / "GSE30019_gene_reox.csv", index=False)
    return pd.DataFrame([met])


def gse240212() -> pd.DataFrame:
    genes, weights, _ = load_theta_panel()
    tar_path = DATA / "GSE240212_RAW.tar"
    dest = DATA / "GSE240212"
    dest.mkdir(exist_ok=True)
    if not list(dest.glob("*matrix.mtx.gz")):
        with tarfile.open(tar_path) as tar:
            tar.extractall(dest)
    rows = []
    for prefix in ["GSM7688224_Tumor-897", "GSM7688225_Tumor-899"]:
        mtx = next(dest.rglob(f"{prefix}_matrix.mtx.gz"))
        feat = next(dest.rglob(f"{prefix}_features.tsv.gz"))
        bc = next(dest.rglob(f"{prefix}_barcodes.tsv.gz"))
        from scipy.io import mmread

        M = mmread(mtx).tocsr().T  # spots x genes
        features = pd.read_csv(feat, sep="\t", header=None)
        symbols = features.iloc[:, 1].astype(str).str.upper().to_numpy() if features.shape[1] > 1 else features.iloc[:, 0].astype(str).str.upper().to_numpy()
        ix = []
        used = []
        w_use = []
        for g, w in zip(genes, weights):
            hits = np.flatnonzero(symbols == g.upper())
            if hits.size:
                ix.append(int(hits[0]))
                used.append(g)
                w_use.append(w)
        X = np.asarray(M[:, ix].todense() if sparse.issparse(M) else M[:, ix], dtype=np.float64)
        logx = log1p_cp10k(X)
        w_use = np.asarray(w_use, dtype=np.float64)
        w_use = w_use / w_use.sum()
        score = module_score(logx, w_use)
        barcodes = pd.read_csv(bc, header=None)[0].astype(str)
        pos = next(dest.rglob(f"{prefix}_tissue_positions_list.csv.gz"))
        spots = pd.read_csv(pos, header=None)
        spots.columns = ["barcode", "in_tissue", "row", "col", "pxl_row", "pxl_col"][: spots.shape[1]]
        spots["barcode"] = spots["barcode"].astype(str)
        df = pd.DataFrame({"barcode": barcodes, "hif_module": score, "tumor": prefix})
        df = df.merge(spots, on="barcode", how="left")
        df.to_csv(OUT / f"{prefix}_spot_scores.csv", index=False)
        in_t = df["in_tissue"].eq(1) if "in_tissue" in df.columns else np.ones(len(df), bool)
        s = df.loc[in_t, "hif_module"].to_numpy()
        # spatial autocorrelation as a sanity check that HIF is spatially structured (not noise)
        rows.append(
            {
                "dataset": "GSE240212",
                "split": f"{prefix} in-tissue spots",
                "n": int(in_t.sum()),
                "median_module": float(np.median(s)),
                "iqr_module": float(np.subtract(*np.percentile(s, [75, 25]))),
                "genes_found": ",".join(used),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    parts = []
    print("=== GSE292771 (priority 5: cancer sc, labeled O2 + HIF KO) ===")
    tab, cells = gse292771()
    cells.to_csv(OUT / "GSE292771_cell_scores.csv", index=False)
    print(tab.to_string(index=False))
    parts.append(tab)
    print("=== GSE30019 (Lai bulk reox time course; gene-set dynamics) ===")
    t19 = gse30019()
    print(t19.to_string(index=False))
    parts.append(t19)
    print("=== GSE240212 (Godet Visium; GFP images not in MTX — HIF spatial structure) ===")
    tvis = gse240212()
    print(tvis.to_string(index=False))
    parts.append(tvis)
    seurat = DATA / "GSE200207_extracted_counts.npz"
    if seurat.exists():
        print("=== GSE200207 extracted ===")
    else:
        print("GSE200207/GSE296547 Seurat RDS need R/Seurat dump (see extract_seurat.R)")
    print("partial tables computed; run robustness.py for the canonical summary.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
