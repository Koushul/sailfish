"""5' GEM-X V(D)J assembly for TCR and BCR.

Chemistry follows Chromium GEM-X Single Cell 5' v3 (CG000733): Read 1 is a
16 bp barcode, a 12 bp UMI, then the 13 bp TSO ``TTTCTTATATGGG``. The barcode
inclusion list is the Single Cell 5' v3 list ``3M-5pgex-jan-2023``.

Contig annotation follows Cell Ranger's V(D)J algorithm: CDR3 runs from the
conserved cysteine through the first residue of an FGXG/WGXG motif (26–80 nt),
and a contig is productive only when it spans a V start codon and the end of J
in frame with no stop. Cell calling uses the three documented filters
(productive contig, at least three UMIs with two read pairs, and the N50
read-per-UMI floor). Clonotypes group cells that share productive CDR3
nucleotide sequences.
"""

from __future__ import annotations

import csv
import gzip
import json
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .quant import ensure_chemistry

TSO = "TTTCTTATATGGG"
BARCODE_LEN = 16
UMI_LEN = 12
CHEMISTRY = "10xv3-5p"
GEOMETRY = "1{b[16]u[12]x:}2{r:}"
WHITELIST_NAME = "3M-5pgex-jan-2023.txt.gz"
WHITELIST_URL = (
    "https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/"
    "REF-10x_barcodes/3M-5pgex-jan-2023.txt.gz"
)
MIN_CDR3_NT = 26
MAX_CDR3_NT = 80
LEFT_MOTIFS = ("LQPEDSAVYY", "VEASQTGTYF", "ATSGQASLYL")
RIGHT_MOTIFS = ("LTFG.GTRVTV", "LIWG.GSKLSI")
TCR_CHAINS = ("TRA", "TRB", "TRG", "TRD")
BCR_CHAINS = ("IGH", "IGK", "IGL")
_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODON = {
    a + b + c: _AA[i * 16 + j * 4 + k]
    for i, a in enumerate("TCAG")
    for j, b in enumerate("TCAG")
    for k, c in enumerate("TCAG")
}
ADAPTER = "AGATCGGAAGAGC"


@dataclass
class Segment:
    gene: str
    region: str
    chain: str
    receptor: str
    seq: str
    anchor: int | None = None


@dataclass
class Reference:
    segments: list[Segment]
    v_by_chain: dict[str, list[Segment]]
    j_by_chain: dict[str, list[Segment]]
    d_by_chain: dict[str, list[Segment]]
    c_by_chain: dict[str, list[Segment]]
    v_kmers: dict[str, list[tuple[str, str]]]
    j_kmers: dict[str, list[tuple[str, str]]]
    d_kmers: dict[str, list[tuple[str, str]]]
    c_kmers: dict[str, list[tuple[str, str]]]


@dataclass
class Cdr3Hit:
    nt: str
    aa: str
    start: int
    end: int
    flank: int


@dataclass
class Group:
    pairs: set[int] = field(default_factory=set)
    umi_pairs: dict[str, set[int]] = field(default_factory=dict)
    seqs: list[tuple[str, int, int]] = field(default_factory=list)


def revcomp(seq: str) -> str:
    return seq.translate(_RC)[::-1].upper()


def translate(seq: str) -> str:
    return "".join(_CODON.get(seq[i : i + 3], "X") for i in range(0, len(seq) - 2, 3))


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def _motif_score(seq: str, motifs: tuple[str, ...]) -> int:
    if not seq or len(seq) != len(motifs[0]):
        return 0
    score = 0
    for i, aa in enumerate(seq):
        if any(m[i] == aa or m[i] == "." for m in motifs):
            score += 1
    return score


def flank_score(prot: str, cys_aa: int, end_aa: int) -> int:
    left = prot[cys_aa - 10 : cys_aa] if cys_aa >= 10 else ""
    right = prot[end_aa + 1 : end_aa + 12]
    return _motif_score(left, LEFT_MOTIFS) + _motif_score(right, RIGHT_MOTIFS)


def find_cdr3(seq: str) -> Cdr3Hit | None:
    """CDR3 from the cysteine through the first residue of FGXG or WGXG."""
    seq = seq.upper()
    best: tuple[int, int, int, Cdr3Hit] | None = None
    for frame in range(3):
        prot: list[str] = []
        pos: list[int] = []
        for i in range(frame, len(seq) - 2, 3):
            prot.append(_CODON.get(seq[i : i + 3], "X"))
            pos.append(i)
        for i in range(len(prot) - 3):
            tetra = prot[i] + prot[i + 1] + prot[i + 2] + prot[i + 3]
            if tetra[0] not in "FW" or tetra[1] != "G" or tetra[3] != "G":
                continue
            c_idx = None
            for c in range(i - 1, -1, -1):
                if prot[c] == "*":
                    break
                if prot[c] != "C":
                    continue
                nt_len = pos[i] + 3 - pos[c]
                if nt_len > MAX_CDR3_NT:
                    break
                if nt_len >= MIN_CDR3_NT:
                    c_idx = c
                    break
            if c_idx is None or c_idx == 0 or prot[c_idx - 1] not in "YFL":
                continue
            aa = "".join(prot[c_idx : i + 1])
            if "*" in aa or not aa.startswith("C"):
                continue
            start = pos[c_idx]
            end = pos[i] + 3
            hit = Cdr3Hit(seq[start:end], aa, start, end, flank_score("".join(prot), c_idx, i))
            rank = (hit.flank, start, len(aa))
            if best is None or rank > best[:3]:
                best = (rank[0], rank[1], rank[2], hit)
    return None if best is None else best[3]


def _region_code(region_type: str) -> str | None:
    return {
        "L-REGION+V-REGION": "V",
        "D-REGION": "D",
        "J-REGION": "J",
        "C-REGION": "C",
        "5'UTR": "U",
    }.get(region_type)


def _v_anchor(seq: str) -> int | None:
    prot = translate(seq)
    best = None
    for i, aa in enumerate(prot):
        if aa != "C" or i < 10:
            continue
        if i < len(prot) - 25:
            continue
        left = prot[i - 10 : i]
        score = _motif_score(left, LEFT_MOTIFS)
        rank = (score, i)
        if best is None or rank > best[0]:
            best = (rank, i * 3)
    if best is None and "C" in prot:
        return prot.rfind("C") * 3
    return None if best is None else best[1]


def _j_anchor(seq: str) -> int | None:
    best = None
    for frame in range(3):
        prot = translate(seq[frame:])
        for i in range(len(prot) - 3):
            tetra = prot[i : i + 4]
            if tetra[0] in "FW" and tetra[1] == "G" and tetra[3] == "G":
                if "*" in prot[:i]:
                    continue
                nt = frame + i * 3
                if best is None or nt < best:
                    best = nt
    return best


def parse_regions_fa(path: Path) -> Reference:
    segments: list[Segment] = []
    header = None
    buf: list[str] = []

    def flush(head: str, seq: str) -> None:
        parts = head[1:].strip().split("|")
        if len(parts) < 6:
            return
        code = _region_code(parts[3])
        if code not in {"V", "D", "J", "C"}:
            return
        chain = parts[5]
        receptor = "TCR" if chain.startswith("TR") else "BCR" if chain.startswith("IG") else ""
        anchor = None
        if code == "V":
            anchor = _v_anchor(seq)
        elif code == "J":
            anchor = _j_anchor(seq)
        segments.append(Segment(parts[2], code, chain, receptor, seq, anchor))

    with Path(path).open() as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    flush(header, "".join(buf))
                header = line
                buf = []
            else:
                buf.append(line.strip().upper())
        if header is not None:
            flush(header, "".join(buf))

    v_by: dict[str, list[Segment]] = defaultdict(list)
    j_by: dict[str, list[Segment]] = defaultdict(list)
    d_by: dict[str, list[Segment]] = defaultdict(list)
    c_by: dict[str, list[Segment]] = defaultdict(list)
    v_kmers: dict[str, list[tuple[str, str]]] = defaultdict(list)
    j_kmers: dict[str, list[tuple[str, str]]] = defaultdict(list)
    d_kmers: dict[str, list[tuple[str, str]]] = defaultdict(list)
    c_kmers: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for seg in segments:
        if seg.region == "V" and seg.anchor is not None and seg.anchor >= 30:
            v_by[seg.chain].append(seg)
            sig = seg.seq[seg.anchor - 60 : seg.anchor]
            for i in range(0, len(sig) - 11):
                v_kmers[sig[i : i + 12]].append((seg.chain, seg.gene))
        elif seg.region == "J" and seg.anchor is not None:
            j_by[seg.chain].append(seg)
            sig = seg.seq[seg.anchor : seg.anchor + 40]
            for i in range(0, max(0, len(sig) - 11)):
                j_kmers[sig[i : i + 12]].append((seg.chain, seg.gene))
        elif seg.region == "D" and len(seg.seq) >= 8:
            d_by[seg.chain].append(seg)
            for i in range(0, len(seg.seq) - 7):
                d_kmers[seg.seq[i : i + 8]].append((seg.chain, seg.gene))
        elif seg.region == "C":
            c_by[seg.chain].append(seg)
            sig = seg.seq[:80]
            for i in range(0, max(0, len(sig) - 11)):
                c_kmers[sig[i : i + 12]].append((seg.chain, seg.gene))
    return Reference(segments, v_by, j_by, d_by, c_by, v_kmers, j_kmers, d_kmers, c_kmers)


def resolve_regions_fa(reference: Path) -> Path:
    reference = Path(reference)
    if reference.is_file():
        return reference
    cand = reference / "fasta" / "regions.fa"
    if cand.exists():
        return cand
    raise FileNotFoundError(f"No fasta/regions.fa under {reference}")


def load_whitelist(path: Path) -> set[str]:
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    out: set[str] = set()
    with opener(path, "rt") as fh:
        for line in fh:
            bc = line.strip().split()[0].upper()
            if bc:
                out.add(bc)
    if not out:
        raise RuntimeError(f"Empty barcode inclusion list: {path}")
    return out


def correct_barcode(bc: str, whitelist: set[str]) -> str | None:
    if bc in whitelist:
        return bc
    if bc.count("N") > 1:
        return None
    hit = None
    for i, base in enumerate(bc):
        for alt in "ACGT":
            if alt == base:
                continue
            nb = bc[:i] + alt + bc[i + 1 :]
            if nb in whitelist:
                if hit is not None and hit != nb:
                    return None
                hit = nb
    return hit


def r1_insert(seq: str) -> tuple[str, str, str]:
    seq = seq.upper()
    bc = seq[:BARCODE_LEN]
    umi = seq[BARCODE_LEN : BARCODE_LEN + UMI_LEN]
    window = seq[28:41]
    if len(window) == len(TSO) and hamming(window, TSO) <= 1:
        insert = seq[41:]
    else:
        at = seq.find(TSO, 20, 55)
        insert = seq[at + len(TSO) :] if at >= 0 else seq[41:]
    return bc, umi, _trim_sense(insert)


def _trim_sense(seq: str) -> str:
    if seq.startswith(TSO):
        seq = seq[len(TSO) :]
    at = seq.find(ADAPTER)
    if at >= 0:
        seq = seq[:at]
    return seq


def _vote(seq: str, index: dict[str, list[tuple[str, str]]], k: int, chain_ok: set[str]) -> Counter:
    votes: Counter = Counter()
    if len(seq) < k:
        return votes
    for i in range(0, len(seq) - k + 1):
        for chain, gene in index.get(seq[i : i + k], ()):
            if chain in chain_ok:
                votes[(chain, gene)] += 1
    return votes


def _gene_seq(ref: Reference, chain: str, gene: str, region: str) -> Segment | None:
    table = {"V": ref.v_by_chain, "J": ref.j_by_chain, "D": ref.d_by_chain, "C": ref.c_by_chain}[region]
    for seg in table.get(chain, ()):
        if seg.gene == gene:
            return seg
    return None


def annotate_read(seq: str, hit: Cdr3Hit, ref: Reference, chains: set[str]) -> dict | None:
    up = seq[max(0, hit.start - 60) : hit.start]
    down = seq[hit.end - 3 : hit.end + 36]
    v_votes = _vote(up, ref.v_kmers, 12, chains)
    j_votes = _vote(down, ref.j_kmers, 12, chains)
    if not v_votes or not j_votes:
        return None
    (v_chain, v_gene), v_n = v_votes.most_common(1)[0]
    j_cands = [(cj, n) for cj, n in j_votes.items() if cj[0] == v_chain]
    if not j_cands:
        return None
    (j_chain, j_gene), j_n = max(j_cands, key=lambda x: x[1])
    if v_n < 2 or j_n < 1:
        return None
    v = _gene_seq(ref, v_chain, v_gene, "V")
    j = _gene_seq(ref, j_chain, j_gene, "J")
    if v is None or j is None or v.anchor is None or j.anchor is None:
        return None
    if hit.start < 30 or v.anchor < 30:
        return None
    if hamming(seq[hit.start - 30 : hit.start], v.seq[v.anchor - 30 : v.anchor]) > 8:
        return None
    mid = seq[hit.start + 3 : hit.end - 3]
    d_votes = _vote(mid, ref.d_kmers, 8, {v_chain})
    d_gene = d_votes.most_common(1)[0][0][1] if d_votes else ""
    tail = seq[hit.end - 3 :]
    c_votes = _vote(tail, ref.c_kmers, 12, {v_chain})
    c_gene = c_votes.most_common(1)[0][0][1] if c_votes else ""
    contig = _build_contig(seq, hit, v, j)
    flags = _productivity(contig, v, j, hit)
    return {
        "chain": v_chain,
        "v_gene": v_gene,
        "d_gene": d_gene,
        "j_gene": j_gene,
        "c_gene": c_gene,
        "cdr3": hit.aa,
        "cdr3_nt": hit.nt,
        "sequence": contig,
        "length": len(contig),
        "full_length": flags["full_length"],
        "productive": flags["productive"],
        "flank": hit.flank,
    }


def _build_contig(seq: str, hit: Cdr3Hit, v: Segment, j: Segment) -> str:
    assert v.anchor is not None and j.anchor is not None
    covered = seq[max(0, hit.start - 80) : hit.start]
    tail = v.seq[: v.anchor]
    if hit.start >= v.anchor:
        start = hit.start - v.anchor
        body = seq[start:]
    elif len(covered) < len(tail):
        body = tail[: len(tail) - len(covered)] + covered + seq[hit.start :]
    else:
        body = tail + seq[hit.start :]
    f_on = len(body) - (len(seq) - (hit.end - 3))
    j_from_f = j.seq[j.anchor :]
    have = len(body) - f_on
    if 0 <= f_on <= len(body) and have < len(j_from_f):
        body = body + j_from_f[have:]
    return body


def _productivity(contig: str, v: Segment, j: Segment, hit: Cdr3Hit) -> dict:
    assert v.anchor is not None and j.anchor is not None
    has_start = contig.startswith("ATG") or (len(contig) >= 3 and contig[:3] == "ATG")
    v_start = 0 if contig.startswith("ATG") else -1
    cdr3_at = contig.find(hit.nt)
    in_frame = has_start and cdr3_at >= 0 and cdr3_at % 3 == 0
    j_end = cdr3_at + len(hit.nt) - 3 + (len(j.seq) - j.anchor) if cdr3_at >= 0 else -1
    j_full = j_end >= cdr3_at and (j.anchor + (j_end - (cdr3_at + len(hit.nt) - 3))) >= len(j.seq) - 3
    no_stop = False
    if in_frame and j_end > v_start:
        no_stop = all(translate(contig[i : i + 3]) != "*" for i in range(0, min(j_end, len(contig) - 2), 3))
    span = j_end - v_start if j_end > 0 else 0
    delta = (len(v.seq) + len(j.seq)) - span
    lo = -55 if v.chain == "IGH" else -25
    structure = lo <= delta <= 25
    order_ok = cdr3_at >= 0
    full_length = has_start and j_full and order_ok
    productive = bool(full_length and in_frame and no_stop and structure and hit.aa.startswith("C") and "*" not in hit.aa)
    return {"full_length": full_length, "productive": productive}


def iter_read_pairs(reads1: list[str], reads2: list[str], max_reads: int | None = None):
    n = 0
    for p1, p2 in zip(reads1, reads2):
        with gzip.open(p1, "rt") as a, gzip.open(p2, "rt") as b:
            while True:
                if not a.readline():
                    break
                s1 = a.readline().strip()
                a.readline()
                a.readline()
                b.readline()
                s2 = b.readline().strip()
                b.readline()
                b.readline()
                if not s1 or not s2:
                    break
                n += 1
                yield s1, s2
                if max_reads is not None and n >= max_reads:
                    return


def resolve_fastqs(lib: dict | None) -> dict | None:
    if not lib:
        return None
    lib = dict(lib)
    r1 = lib.get("reads1") or []
    r2 = lib.get("reads2") or []
    if isinstance(r1, str):
        r1 = [r1]
    if isinstance(r2, str):
        r2 = [r2]
    if r1 and r2:
        lib["reads1"], lib["reads2"] = list(r1), list(r2)
        return lib
    folder = lib.get("fastqs")
    if not folder or not Path(folder).is_dir():
        return lib
    r1 = sorted(str(p) for p in Path(folder).glob("*_R1_001.fastq.gz"))
    r2 = sorted(str(p) for p in Path(folder).glob("*_R2_001.fastq.gz"))
    lib["reads1"], lib["reads2"] = r1, r2
    return lib


def _keep_seq(group: Group, seq: str, start: int, flank: int) -> None:
    group.seqs.append((seq, start, flank))
    if len(group.seqs) > 6:
        group.seqs.sort(key=lambda x: (x[2], len(x[0])), reverse=True)
        del group.seqs[6:]


def count_barcodes(reads1, reads2, whitelist: set[str], max_reads: int | None) -> Counter:
    counts: Counter = Counter()
    n = 0
    for s1, _s2 in iter_read_pairs(reads1, reads2, max_reads):
        n += 1
        bc, _umi, _ins = r1_insert(s1)
        if "N" in bc:
            fixed = correct_barcode(bc, whitelist)
        else:
            fixed = bc if bc in whitelist else correct_barcode(bc, whitelist)
        if fixed:
            counts[fixed] += 1
        if n % 5000000 == 0:
            print(f"vdj barcode pass {n:,} reads, {len(counts):,} barcodes", flush=True)
    print(f"vdj barcode pass done {n:,} reads, {len(counts):,} barcodes", flush=True)
    return counts


def collect_groups(reads1, reads2, whitelist, keep: set[str], ref: Reference, chains: set[str], max_reads: int | None):
    groups: dict[tuple[str, str], Group] = {}
    n = 0
    n_hit = 0
    for s1, s2 in iter_read_pairs(reads1, reads2, max_reads):
        n += 1
        bc, umi, insert = r1_insert(s1)
        fixed = bc if bc in whitelist else correct_barcode(bc, whitelist)
        if fixed not in keep or "N" in umi:
            continue
        seqs = []
        if len(insert) >= 60:
            seqs.append(insert)
        r2 = _trim_sense(revcomp(s2))
        if len(r2) >= 60:
            seqs.append(r2)
        seen = set()
        for seq in seqs:
            hit = find_cdr3(seq)
            if hit is None or hit.nt in seen:
                continue
            ann = annotate_read(seq, hit, ref, chains)
            if ann is None:
                continue
            seen.add(hit.nt)
            n_hit += 1
            key = (fixed, hit.nt)
            group = groups.get(key)
            if group is None:
                group = Group()
                groups[key] = group
            group.pairs.add(n)
            group.umi_pairs.setdefault(umi, set()).add(n)
            _keep_seq(group, seq, hit.start, hit.flank)
        if n % 2000000 == 0:
            print(f"vdj annotate pass {n:,} reads, {len(groups):,} CDR3 groups, {n_hit:,} hits", flush=True)
    print(f"vdj annotate pass done {n:,} reads, {len(groups):,} groups", flush=True)
    return groups, n


def _collapse_umis(umi_pairs: dict[str, set[int]]) -> dict[str, set[int]]:
    if len(umi_pairs) < 2:
        return {u: set(ps) for u, ps in umi_pairs.items()}
    parent = {u: u for u in umi_pairs}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def size(x: str) -> int:
        return len(umi_pairs[find(x)])

    for u in list(umi_pairs):
        for i, base in enumerate(u):
            for alt in "ACGT":
                if alt == base:
                    continue
                nb = u[:i] + alt + u[i + 1 :]
                if nb not in parent or find(nb) == find(u):
                    continue
                if size(nb) >= size(u):
                    parent[find(u)] = find(nb)
    merged: dict[str, set[int]] = {}
    for u, pairs in umi_pairs.items():
        merged.setdefault(find(u), set()).update(pairs)
    return merged


def contigs_from_groups(groups: dict[tuple[str, str], Group], ref: Reference, chains: set[str], receptor: str) -> list[dict]:
    rows = []
    for (bc, _nt), group in groups.items():
        if not group.seqs:
            continue
        seq, _start, _flank = max(group.seqs, key=lambda x: (x[2], len(x[0])))
        hit = find_cdr3(seq)
        if hit is None:
            continue
        ann = annotate_read(seq, hit, ref, chains)
        if ann is None or ann["chain"] not in chains:
            continue
        merged = _collapse_umis(group.umi_pairs)
        depths = [len(ps) for ps in merged.values()]
        rows.append(
            {
                **ann,
                "barcode": bc,
                "receptor": receptor,
                "reads": len(group.pairs),
                "umis": len(merged),
                "filtered_umis": sum(1 for d in depths if d >= 2),
                "max_umi_reads": max(depths) if depths else 0,
            }
        )
    return rows


def _n50(values: list[int]) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    total = sum(ordered)
    acc = 0
    for v in ordered:
        acc += v
        if acc >= total / 2:
            return float(v)
    return float(ordered[-1])


def call_cells(rows: list[dict]) -> None:
    by_bc: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        by_bc[row["barcode"]].append(row)
    umi_depths = []
    for group in by_bc.values():
        for row in group:
            if row["filtered_umis"]:
                umi_depths.extend([row["max_umi_reads"]] * row["filtered_umis"])
    n50 = _n50(umi_depths) if umi_depths else 0.0
    floor = 0.03 * n50
    for bc, group in by_bc.items():
        productive = [r for r in group if r["productive"]]
        drop = False
        by_chain: dict[str, list[dict]] = defaultdict(list)
        for r in productive:
            by_chain[r["chain"]].append(r)
        for chain_rows in by_chain.values():
            strong = [r for r in chain_rows if r["umis"] >= 2]
            cdrs = {r["cdr3_nt"] for r in strong}
            vs = {r["v_gene"] for r in strong}
            if len(cdrs) > 1 and len(vs) == 1:
                drop = True
        n_filt = sum(r["filtered_umis"] for r in group)
        max_depth = max((r["max_umi_reads"] for r in group), default=0)
        only = len(productive) == 1
        junction_ok = (not only) or (productive and productive[0]["umis"] > 1)
        is_cell = (
            bool(productive)
            and junction_ok
            and n_filt >= 3
            and max_depth >= floor
            and not drop
            and len(productive) <= 4
        )
        for r in group:
            r["is_cell"] = is_cell
            r["high_confidence"] = bool(r["productive"] and r["umis"] >= 2 and not drop)
            r["raw_clonotype_id"] = ""
            r["raw_consensus_id"] = ""


def assign_clonotypes(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    cells: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        if row.get("is_cell") and row.get("productive"):
            cells[row["barcode"]].append(row)
    groups: dict[tuple, list[str]] = defaultdict(list)
    for bc, group in cells.items():
        key = tuple(sorted((r["chain"], r["cdr3_nt"]) for r in group))
        if key:
            groups[key].append(bc)
    ordered = sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    clonotypes = []
    bc_to_id = {}
    for i, (key, bcs) in enumerate(ordered, start=1):
        cid = f"clonotype{i}"
        for bc in bcs:
            bc_to_id[bc] = cid
        example = cells[bcs[0]]
        aa = ";".join(f"{r['chain']}:{r['cdr3']}" for r in sorted(example, key=lambda r: r["chain"]))
        nt = ";".join(f"{r['chain']}:{r['cdr3_nt']}" for r in sorted(example, key=lambda r: r["chain"]))
        clonotypes.append(
            {
                "clonotype_id": cid,
                "frequency": len(bcs),
                "proportion": 0.0,
                "cdr3s_aa": aa,
                "cdr3s_nt": nt,
                "key": key,
            }
        )
    n_cells = sum(c["frequency"] for c in clonotypes) or 1
    for c in clonotypes:
        c["proportion"] = c["frequency"] / n_cells
    for row in rows:
        cid = bc_to_id.get(row["barcode"], "")
        if row.get("is_cell") and row.get("productive"):
            row["raw_clonotype_id"] = cid
    consensus = []
    for c in clonotypes:
        members = []
        for bc in groups[c["key"]]:
            members.extend(cells[bc])
        by_chain: dict[str, list[dict]] = defaultdict(list)
        for r in members:
            by_chain[r["chain"]].append(r)
        for chain, chain_rows in sorted(by_chain.items()):
            top = max(chain_rows, key=lambda r: (r["umis"], r["reads"]))
            cons_id = f"{c['clonotype_id']}_consensus_{chain}"
            for r in chain_rows:
                r["raw_consensus_id"] = cons_id
            consensus.append(
                {
                    "clonotype_id": c["clonotype_id"],
                    "consensus_id": cons_id,
                    "length": top["length"],
                    "chain": chain,
                    "v_gene": top["v_gene"],
                    "d_gene": top["d_gene"],
                    "j_gene": top["j_gene"],
                    "c_gene": top["c_gene"],
                    "full_length": top["full_length"],
                    "productive": top["productive"],
                    "cdr3": top["cdr3"],
                    "cdr3_nt": top["cdr3_nt"],
                    "sequence": top["sequence"],
                    "reads": sum(r["reads"] for r in chain_rows),
                    "umis": sum(r["umis"] for r in chain_rows),
                }
            )
    return clonotypes, consensus


def _write_csv(path: Path, rows: list[dict], columns: list[str], delimiter: str = ",") -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore", delimiter=delimiter)
        w.writeheader()
        for row in rows:
            w.writerow(row)


CONTIG_COLS = [
    "barcode",
    "is_cell",
    "contig_id",
    "high_confidence",
    "length",
    "chain",
    "v_gene",
    "d_gene",
    "j_gene",
    "c_gene",
    "full_length",
    "productive",
    "cdr3",
    "cdr3_nt",
    "reads",
    "umis",
    "raw_clonotype_id",
    "raw_consensus_id",
]


def write_library(outdir: Path, receptor: str, rows: list[dict], n_reads: int, n_barcodes_kept: int) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    rows.sort(key=lambda r: (r["barcode"], r["chain"], -r["umis"]))
    for i, row in enumerate(rows, start=1):
        row["contig_id"] = f"{row['barcode']}_contig_{i}"
    clonotypes, consensus = assign_clonotypes(rows)
    filtered = [r for r in rows if r.get("is_cell") and r.get("productive")]
    _write_csv(outdir / "all_contig_annotations.csv", rows, CONTIG_COLS)
    _write_csv(outdir / "filtered_contig_annotations.csv", filtered, CONTIG_COLS)
    _write_csv(
        outdir / "clonotypes.csv",
        clonotypes,
        ["clonotype_id", "frequency", "proportion", "cdr3s_aa", "cdr3s_nt"],
    )
    _write_csv(
        outdir / "consensus_annotations.csv",
        consensus,
        [
            "clonotype_id",
            "consensus_id",
            "length",
            "chain",
            "v_gene",
            "d_gene",
            "j_gene",
            "c_gene",
            "full_length",
            "productive",
            "cdr3",
            "cdr3_nt",
            "reads",
            "umis",
        ],
    )
    airr = []
    for row in filtered:
        airr.append(
            {
                "sequence_id": row["contig_id"],
                "sequence": row["sequence"],
                "productive": str(row["productive"]).upper(),
                "locus": row["chain"],
                "v_call": row["v_gene"],
                "d_call": row["d_gene"],
                "j_call": row["j_gene"],
                "c_call": row["c_gene"],
                "junction": row["cdr3_nt"],
                "junction_aa": row["cdr3"],
                "umi_count": row["umis"],
                "consensus_count": row["reads"],
                "cell_id": row["barcode"],
                "clone_id": row["raw_clonotype_id"],
            }
        )
    _write_csv(
        outdir / "airr_rearrangement.tsv",
        airr,
        [
            "sequence_id",
            "sequence",
            "productive",
            "locus",
            "v_call",
            "d_call",
            "j_call",
            "c_call",
            "junction",
            "junction_aa",
            "umi_count",
            "consensus_count",
            "cell_id",
            "clone_id",
        ],
        delimiter="\t",
    )
    with (outdir / "filtered_contig.fasta").open("w") as fh:
        for row in filtered:
            fh.write(f">{row['contig_id']}\n{row['sequence']}\n")
    cells = {r["barcode"] for r in filtered}
    paired = 0
    for bc in cells:
        chains = {r["chain"] for r in filtered if r["barcode"] == bc}
        if receptor == "TCR" and "TRA" in chains and "TRB" in chains:
            paired += 1
        if receptor == "BCR" and "IGH" in chains and ({"IGK", "IGL"} & chains):
            paired += 1
    metrics = {
        "receptor": receptor,
        "number_of_reads": n_reads,
        "barcodes_assembled": n_barcodes_kept,
        "number_of_cells": len(cells),
        "number_of_paired_cells": paired,
        "number_of_clonotypes": len(clonotypes),
        "number_of_productive_contigs": len(filtered),
        "chemistry": CHEMISTRY,
        "barcode_whitelist": WHITELIST_NAME,
    }
    with (outdir / "metrics_summary.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["metric", "value"])
        for k, v in metrics.items():
            w.writerow([k, v])
    (outdir / "vdj_summary.json").write_text(json.dumps(metrics, indent=2) + "\n")
    return metrics


def run_vdj_library(lib: dict, receptor: str, ref: Reference, whitelist: set[str], outdir: Path, vdj_cfg: dict) -> dict:
    chains = set(TCR_CHAINS if receptor == "TCR" else BCR_CHAINS)
    max_reads = vdj_cfg.get("max_reads")
    min_pairs = int(lib.get("min_read_pairs", vdj_cfg.get("min_read_pairs", 6)))
    counts = count_barcodes(lib["reads1"], lib["reads2"], whitelist, max_reads)
    keep = {bc for bc, n in counts.items() if n >= min_pairs}
    print(f"vdj {receptor}: {len(keep):,} barcodes with ≥{min_pairs} read pairs", flush=True)
    groups, n_reads = collect_groups(lib["reads1"], lib["reads2"], whitelist, keep, ref, chains, max_reads)
    rows = contigs_from_groups(groups, ref, chains, receptor)
    call_cells(rows)
    return write_library(outdir, receptor, rows, n_reads, len(keep))


def _register_chemistry(whitelist_txt: Path) -> None:
    try:
        ensure_chemistry(CHEMISTRY, GEOMETRY, "fw", whitelist_txt)
    except Exception as exc:
        print(f"vdj: chemistry {CHEMISTRY} not registered ({exc})", flush=True)


def ensure_whitelist(path: Path | None, cache: Path) -> Path:
    if path and Path(path).exists():
        return Path(path)
    dest = cache / WHITELIST_NAME
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    cache.mkdir(parents=True, exist_ok=True)
    import urllib.request

    print(f"vdj: downloading {WHITELIST_NAME}", flush=True)
    urllib.request.urlretrieve(WHITELIST_URL, dest)
    return dest


def run_vdj(cfg: dict, outdir: Path) -> dict | None:
    vdj = cfg.get("vdj")
    if not vdj:
        return None
    tcr = resolve_fastqs(vdj.get("tcr"))
    bcr = resolve_fastqs(vdj.get("bcr"))

    def ready(lib):
        return bool(lib and lib.get("reads1") and lib.get("reads2") and all(Path(p).exists() for p in lib["reads1"] + lib["reads2"]))

    if not ready(tcr) and not ready(bcr):
        print("vdj: no TCR/BCR FASTQs found; skipping")
        return {"skipped": True, "reason": "FASTQs missing"}
    reference = vdj.get("reference")
    if not reference:
        raise ValueError("vdj.reference is required (10x refdata-cellranger-vdj directory)")
    ref = parse_regions_fa(resolve_regions_fa(Path(reference)))
    wl_path = ensure_whitelist(Path(vdj["whitelist"]) if vdj.get("whitelist") else None, outdir / "vdj")
    whitelist = load_whitelist(wl_path)
    plain = wl_path
    if wl_path.suffix == ".gz":
        plain = outdir / "vdj" / "3M-5pgex-jan-2023.txt"
        if not plain.exists():
            with gzip.open(wl_path, "rt") as src, plain.open("w") as dst:
                dst.write(src.read())
    _register_chemistry(plain)
    out: dict = {"chemistry": CHEMISTRY, "whitelist": str(wl_path), "reference": str(reference), "libraries": {}}
    if ready(tcr):
        out["libraries"]["TCR"] = run_vdj_library(tcr, "TCR", ref, whitelist, outdir / "vdj" / "tcr", vdj)
    if ready(bcr):
        out["libraries"]["BCR"] = run_vdj_library(bcr, "BCR", ref, whitelist, outdir / "vdj" / "bcr", vdj)
    (outdir / "vdj" / "vdj_run.json").write_text(json.dumps(out, indent=2) + "\n")
    return out
