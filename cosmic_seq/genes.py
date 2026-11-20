"""Full-length genes around the amplicon products of a genome.

A gene is found from the best product of one of its anchor amplicons (step 03), not from annotation names: the GFF3
feature of the gene's type that overlaps the product is exported in full, in gene orientation. Without an overlapping
feature (no GFF delivered, pseudogene, unannotated assembly) the H37Rv gene is aligned to the contig around the product
and the record is cut at the alignment ends (`source=homology`). A record always contains the whole amplicon; the gene
lies at `gene_start..gene_end` (1-based) inside it.

Also the pairwise alignment used here and in step 05 (Biopython PairwiseAligner, blastn-like scores).
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Bio.Align import PairwiseAligner

from .insilico import read_fasta, revcomp
from .tables import Record

FEATURE_TYPES = ("CDS", "rRNA")
MIN_HOMOLOGY_COVERAGE = 0.8  # fraction of the H37Rv gene a homology hit must cover
HOMOLOGY_MARGIN = 200        # bp of contig searched beyond the amplicon + H37Rv gene length on each side
GeneRecord = dict[str, "str | int | float"]


@dataclass(frozen=True)
class GeneSpec:
    gene: str
    anchors: tuple[str, ...]  # amplicons of amplicons.tsv, tried in order
    feature_type: str         # CDS / rRNA
    name_pattern: str         # regex on the GFF attributes; a sanity flag only, never used to pick the feature
    h37rv_locus: str          # H37Rv locus tag: homology fallback and identity; empty if the gene is absent from H37Rv
    paralog_locus: str        # H37Rv locus tag of a paralog the anchor primers may also hit; empty if none


def load_genes(path: Path, amplicon_names: set[str]) -> dict[str, GeneSpec]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    specs: dict[str, GeneSpec] = {}
    for row in rows:
        spec = GeneSpec(row["gene"], tuple(a.strip() for a in row["anchor_amplicons"].split(",") if a.strip()),
                        row["feature_type"], row["name_pattern"], row["h37rv_locus"], row["paralog_locus"])
        if spec.gene in specs:
            raise ValueError(f"{path.name}: gene {spec.gene!r} listed twice")
        if spec.feature_type not in FEATURE_TYPES:
            raise ValueError(f"{path.name}: {spec.gene}: feature_type must be one of {FEATURE_TYPES}")
        unknown = [a for a in spec.anchors if a not in amplicon_names]
        if not spec.anchors or unknown:
            raise ValueError(f"{path.name}: {spec.gene}: anchor amplicons not in amplicons.tsv: {unknown or 'none given'}")
        re.compile(spec.name_pattern)
        specs[spec.gene] = spec
    return specs


# --------------------------------------------------------------------------------------------------------- GFF3
@dataclass(frozen=True)
class Feature:
    type: str   # GFF column 3
    contig: str
    start: int  # 1-based, inclusive
    end: int
    strand: str
    attributes: str


def read_gff(path: Path, types: set[str]) -> list[Feature]:
    """Features of the given types; rows sharing an ID (a CDS split by a frameshift) are merged to one span."""
    spans: dict[tuple[str, str], list[str | int]] = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in types:
                continue
            match = re.search(r"(?:^|;)ID=([^;]+)", fields[8])
            key = (fields[2], match.group(1) if match else f"{fields[0]}:{fields[3]}")
            start, end = int(fields[3]), int(fields[4])
            if key in spans:
                span = spans[key]
                span[1], span[2] = min(int(span[1]), start), max(int(span[2]), end)
            else:
                spans[key] = [fields[0], start, end, fields[6], fields[8]]
    return [Feature(t, str(c), int(s), int(e), str(st), str(a)) for (t, _), (c, s, e, st, a) in spans.items()]


def locus_sequence(fasta: Path, gff: Path, locus: str) -> str:
    """Sequence of the gene with this locus tag, 5'->3' (the `gene` feature)."""
    pattern = re.compile(rf"(?:^|;)locus_tag={re.escape(locus)}(?:;|$)")
    feature = next((f for f in read_gff(gff, {"gene"}) if pattern.search(f.attributes)), None)
    if feature is None:
        raise ValueError(f"locus {locus} not in {gff}")
    seq = read_fasta(fasta)[feature.contig][feature.start - 1:feature.end].upper()
    return revcomp(seq) if feature.strand == "-" else seq


# ---------------------------------------------------------------------------------------------------- alignment
def aligner(local: bool = False) -> Any:
    """Blastn-like scores. Global mode has free end gaps: records of different extent are compared on their overlap."""
    a: Any = PairwiseAligner()  # type: ignore[no-untyped-call]
    a.match_score, a.mismatch_score = 2, -3
    a.open_gap_score, a.extend_gap_score = -5, -2
    a.mode = "local" if local else "global"
    if not local:
        a.end_gap_score = 0
    return a


@dataclass(frozen=True)
class Mapping:
    """Positions of sequence A on sequence B (both 0-based) after one alignment."""

    a_to_b: list[int | None]
    score: float
    identity: float     # matches / alignment columns between the first and last aligned position
    reverse: bool       # B was reverse-complemented to align (only when both orientations were tried)


def align(a: str, b: str, both_strands: bool = False, local: bool = False) -> Mapping:
    """Align B to A; with both_strands the better orientation of B is used (positions then refer to revcomp(B))."""
    al = aligner(local)
    candidates = [(False, b)] + ([(True, revcomp(b))] if both_strands else [])
    scored = [(float(al.score(a, seq)), reverse, seq) for reverse, seq in candidates]
    score, reverse, seq = max(scored, key=lambda t: t[0])
    alignment = al.align(a, seq)[0]
    blocks_a, blocks_b = alignment.aligned
    a_to_b: list[int | None] = [None] * len(a)
    matches = aligned = 0
    for (a0, a1), (b0, b1) in zip(blocks_a, blocks_b):
        for i in range(int(a1) - int(a0)):
            a_to_b[int(a0) + i] = int(b0) + i
            matches += a[int(a0) + i] == seq[int(b0) + i]
        aligned += int(a1) - int(a0)
    if aligned:
        span_a = int(blocks_a[-1][1]) - int(blocks_a[0][0])
        span_b = int(blocks_b[-1][1]) - int(blocks_b[0][0])
        columns = aligned + (span_a - aligned) + (span_b - aligned)
    else:
        columns = 1
    return Mapping(a_to_b, score, matches / columns, reverse)


def identity(a: str, b: str) -> float:
    return align(a, b).identity


# ----------------------------------------------------------------------------------------------- gene extraction
@dataclass(frozen=True)
class References:
    """H37Rv sequences of a gene: homology fallback and identity; paralog for the paralog guard."""

    gene: str | None
    paralog: str | None


def _occurrences(contigs: dict[str, str], amplicon: str) -> list[tuple[str, int, int]]:
    """Every exact location (contig, 1-based start, end) of an amplicon sequence, on either strand."""
    found: list[tuple[str, int, int]] = []
    for target in {amplicon, revcomp(amplicon)}:
        for contig, seq in contigs.items():
            pos = seq.find(target)
            while pos != -1:
                found.append((contig, pos + 1, pos + len(target)))
                pos = seq.find(target, pos + 1)
    return sorted(set(found))


def _homology(contig_seq: str, amp_start: int, amp_end: int, ref: str) -> tuple[int, int, str] | None:
    """Gene span (1-based, strand) found by aligning the H37Rv gene to the contig around the amplicon."""
    w0 = max(1, amp_start - len(ref) - HOMOLOGY_MARGIN)
    w1 = min(len(contig_seq), amp_end + len(ref) + HOMOLOGY_MARGIN)
    window = contig_seq[w0 - 1:w1]
    al = aligner(local=True)
    best: tuple[float, bool, str] | None = None
    for reverse, seq in ((False, window), (True, revcomp(window))):
        score = float(al.score(seq, ref))
        if best is None or score > best[0]:
            best = (score, reverse, seq)
    assert best is not None
    _, reverse, seq = best
    alignment = al.align(seq, ref)[0]
    blocks_w, blocks_r = alignment.aligned
    if not len(blocks_w):
        return None
    r0, r1 = int(blocks_r[0][0]), int(blocks_r[-1][1])
    if (r1 - r0) < MIN_HOMOLOGY_COVERAGE * len(ref):
        return None
    # extend by the unaligned ends of the H37Rv gene (diverged first / last codons), within the window
    s0 = max(0, int(blocks_w[0][0]) - r0)
    s1 = min(len(seq), int(blocks_w[-1][1]) + (len(ref) - r1))
    if reverse:  # coordinates on revcomp(window) -> on the plus strand
        s0, s1 = len(seq) - s1, len(seq) - s0
    start, end = w0 + s0, w0 + s1 - 1
    if end < amp_start or start > amp_end:  # hit elsewhere in the window, not at the amplicon
        return None
    return start, end, "-" if reverse else "+"


def extract(spec: GeneSpec, anchor: Record, contigs: dict[str, str], features: list[Feature] | None,
            refs: References) -> list[GeneRecord]:
    """Gene records at every copy of one anchor product (copies merged by step 03 are located again here)."""
    records: list[GeneRecord] = []
    locations = _occurrences(contigs, str(anchor["sequence"])) or [
        (str(anchor["contig"]), int(anchor["start"]), int(anchor["end"]))]
    for contig, amp_start, amp_end in locations:
        contig_seq = contigs[contig]
        overlapping = [f for f in features or [] if f.contig == contig and f.start <= amp_end and f.end >= amp_start]
        name_ok = ""
        if overlapping:
            feature = max(overlapping, key=lambda f: min(f.end, amp_end) - max(f.start, amp_start))
            gene_start, gene_end, strand, source = feature.start, feature.end, feature.strand, "gff"
            name_ok = "t" if re.search(spec.name_pattern, feature.attributes, re.IGNORECASE) else "f"
        elif refs.gene:
            hit = _homology(contig_seq, amp_start, amp_end, refs.gene)
            if hit is None:
                continue
            gene_start, gene_end, strand = hit
            source = "homology"
        else:
            continue
        rec_start, rec_end = min(gene_start, amp_start), max(gene_end, amp_end)
        sequence = contig_seq[rec_start - 1:rec_end].upper()
        if strand == "-":
            sequence = revcomp(sequence)
            offset = rec_end - gene_end + 1
        else:
            offset = gene_start - rec_start + 1
        records.append({"gene": spec.gene, "anchor": str(anchor["amplicon"]), "contig": contig, "start": rec_start,
                        "end": rec_end, "strand": strand, "length": len(sequence), "gene_start": offset,
                        "gene_end": offset + gene_end - gene_start, "source": source, "name_ok": name_ok,
                        "sequence": sequence})
    return records


def genome_genes(fasta: Path, gff: Path | None, specs: dict[str, GeneSpec], anchors: list[Record],
                 refs: dict[str, References]) -> list[GeneRecord]:
    """Distinct gene records of one genome for every spec; a spec whose anchor gave nothing usable gets a
    `{"gene": g, "failed": anchor}` record when an anchor product existed."""
    by_amplicon: dict[str, list[Record]] = {}
    for product in anchors:
        by_amplicon.setdefault(str(product["amplicon"]), []).append(product)
    needed = [s for s in specs.values() if any(a in by_amplicon for a in s.anchors)]
    if not needed:
        return []
    contigs = read_fasta(fasta)
    types = {s.feature_type for s in needed}
    all_features = read_gff(gff, types) if gff else None
    out: list[GeneRecord] = []
    for spec in needed:
        features = [f for f in all_features if f.type == spec.feature_type] if all_features is not None else None
        anchor_name = next(a for a in spec.anchors if a in by_amplicon)
        found: list[GeneRecord] = []
        for anchor in by_amplicon[anchor_name]:
            found += extract(spec, anchor, contigs, features, refs[spec.gene])
        merged: dict[str, GeneRecord] = {}
        for record in found:
            key = str(record["sequence"])
            if key in merged:
                if (merged[key]["contig"], merged[key]["start"]) != (record["contig"], record["start"]):
                    merged[key]["copies"] = int(merged[key]["copies"]) + 1
                continue
            record["copies"] = 1
            merged[key] = record
        if not merged:
            out.append({"gene": spec.gene, "failed": anchor_name})
            continue
        for record in merged.values():
            gene_seq = str(record["sequence"])[int(record["gene_start"]) - 1:int(record["gene_end"])]
            ref = refs[spec.gene]
            ident_ref = identity(ref.gene, gene_seq) if ref.gene else None
            ident_par = identity(ref.paralog, gene_seq) if ref.paralog else None
            record["identity_h37rv"] = round(ident_ref, 4) if ident_ref is not None else ""
            record["identity_paralog"] = round(ident_par, 4) if ident_par is not None else ""
            record["paralog_flag"] = "t" if ident_ref is not None and ident_par is not None and ident_par > ident_ref else "f"
            out.append(record)
    return out

