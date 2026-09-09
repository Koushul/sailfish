from __future__ import annotations

import gzip
from pathlib import Path


def subsample_paired_fastq(
    r1: Path,
    r2: Path,
    out_r1: Path,
    out_r2: Path,
    n_reads: int,
) -> int:
    out_r1.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with gzip.open(r1, "rt") as a, gzip.open(r2, "rt") as b, gzip.open(out_r1, "wt") as oa, gzip.open(out_r2, "wt") as ob:
        while n < n_reads:
            rec1 = [a.readline() for _ in range(4)]
            rec2 = [b.readline() for _ in range(4)]
            if not rec1[0] or not rec2[0]:
                break
            oa.write("".join(rec1))
            ob.write("".join(rec2))
            n += 1
    return n
