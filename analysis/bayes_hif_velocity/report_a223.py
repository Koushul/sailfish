#!/usr/bin/env python3
"""A223 cannot identify a never-hypoxic control. Report counts only; do not fit."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

TUMOR = "/ix1/ylee/kor11/A223/tumor_kinetics/tumor_qc.csv"
NEU = "/ix1/ylee/kor11/A223/neutrophil_hypoxia/neutrophils.csv"
OUT = Path(__file__).resolve().parent / "results" / "a223_dn_counts.md"
MIN_CONTROL = 40


def _gate_table(df: pd.DataFrame, lineage: str) -> pd.DataFrame:
    rows = []
    for (lane, gate), sub in df.groupby(["lane", "sample_id"]):
        rows.append(
            {
                "lineage": lineage,
                "lane": lane,
                "ocm_gate": gate,
                "n": int(len(sub)),
                "median_spliced_umi": float(sub["spliced_umi"].median()) if "spliced_umi" in sub else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    tumor = pd.read_csv(TUMOR)
    neu = pd.read_csv(NEU)
    ttab = _gate_table(tumor, "Tumor")
    ntab = _gate_table(neu, "neutrophil")
    tab = pd.concat([ttab, ntab], ignore_index=True)
    dn = tab[tab["ocm_gate"] == "DN"]
    lines = [
        "# A223 DN counts",
        "",
        "`hypoxia_plus` is Image-iT+ in the GFP channel, the same probe as E15S (not transgenic GFP, not ROS/DCF). `DN` is Image-iT−. The spliced factor needs a control median/MAD; lag pins mean \(\\xi\) on the control. Counts below the minimum (n_control < 40 **per chemistry**) mean persist/reverted is not identified.",
        "",
        "Do not transfer the E14 control location onto this tumor. Do not pool E27 (3′) with E29 (5′) to inflate DN n. Do not write into A223 h5ads.",
        "",
        "## OCM gate counts (QC tables)",
        "",
    ]
    try:
        body = tab.to_markdown(index=False)
    except Exception:
        body = "```\n" + tab.to_string(index=False) + "\n```"
    lines.append(body)
    lines += [
        "",
        "## Skip decision",
        "",
    ]
    for _, r in dn.iterrows():
        ok = r["n"] >= MIN_CONTROL
        lines.append(
            f"- {r['lineage']} {r['lane']} DN n={int(r['n'])} median UMI={r['median_spliced_umi']:.0f}: "
            + ("would fit" if ok else f"skip (need ≥{MIN_CONTROL})")
        )
    lines += [
        "",
        "Pooled DN tumors across lanes = "
        + str(int(dn.loc[dn["lineage"] == "Tumor", "n"].sum()))
        + "; pooled DN neutrophils = "
        + str(int(dn.loc[dn["lineage"] == "neutrophil", "n"].sum()))
        + ". Still below a credible control for MAD scaling even if chemistry is ignored.",
        "",
        "What would be needed: a never-hypoxic (or clearly HIF-off) arm with tens of cells **per 3′/5′ library**, or a different model that does not center \(h\) and \(\\xi\) on DN.",
        "",
    ]
    OUT.write_text("\n".join(lines) + "\n")
    tab.to_csv(OUT.with_name("a223_gate_counts.tsv"), sep="\t", index=False)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
