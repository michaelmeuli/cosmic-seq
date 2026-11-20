"""The two curated tables of the repo: `primers.tsv` (one row per primer) and `amplicons.tsv` (one row per primer pair),
the products of a genome and the per-genome cache of step 03."""
from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from .insilico import Product, Sites, clean_primer, encode, primer_sites, products, read_fasta

DIRECTIONS = ("forward", "reverse")
Record = dict[str, "str | int"]  # one distinct product sequence of a genome (a line of the step-03 cache)


@dataclass(frozen=True)
class Primer:
    name: str
    target: str
    scope: str
    direction: str
    sequence: str  # upper-case IUPAC, 5'->3'
    role: str
    protocol: str
    source_file: str  # relative to HANDBOOK_DIR


@dataclass(frozen=True)
class Amplicon:
    name: str
    target: str
    scope: str
    forward: str
    reverse: str
    documented_bp: int | None  # product length written in the protocol
    use: str
    sequenced_with: str
    protocol: str


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_primers(path: Path) -> dict[str, Primer]:
    primers: dict[str, Primer] = {}
    for row in _rows(path):
        primer = Primer(row["primer"], row["target"], row["scope"], row["direction"], clean_primer(row["sequence"]),
                        row["role"], row["protocol"], row["source_file"])
        if primer.name in primers:
            raise ValueError(f"{path.name}: primer {primer.name!r} listed twice")
        if primer.direction not in DIRECTIONS:
            raise ValueError(f"{path.name}: primer {primer.name!r} has direction {primer.direction!r}")
        primers[primer.name] = primer
    return primers


def load_amplicons(path: Path, primers: dict[str, Primer]) -> dict[str, Amplicon]:
    amplicons: dict[str, Amplicon] = {}
    for row in _rows(path):
        amplicon = Amplicon(row["amplicon"], row["target"], row["scope"], row["forward"], row["reverse"],
                            int(row["documented_bp"]) if row["documented_bp"] else None, row["use"],
                            row["sequenced_with"], row["protocol"])
        if amplicon.name in amplicons:
            raise ValueError(f"{path.name}: amplicon {amplicon.name!r} listed twice")
        for name, direction in ((amplicon.forward, "forward"), (amplicon.reverse, "reverse")):
            if name not in primers:
                raise ValueError(f"{path.name}: amplicon {amplicon.name!r} names unknown primer {name!r}")
            if primers[name].direction != direction:
                raise ValueError(f"{path.name}: amplicon {amplicon.name!r} uses {name!r} as {direction} primer, "
                                 f"primers.tsv says {primers[name].direction}")
        amplicons[amplicon.name] = amplicon
    return amplicons


def genome_fasta(genome_dir: Path) -> Path | None:
    """The genome FASTA of an NCBI datasets directory (<accession>/<accession>_<assembly>_genomic.fna)."""
    return next(iter(sorted(genome_dir.glob("*_genomic.fna"))), None)


def best(found: list[Product]) -> list[Product]:
    """Products with the fewest primer mismatches in total.

    Weaker products of the same pair are other loci (the groEL1 paralog for the hsp65 primers, unrelated genes
    hit with 3 + 3 mismatches). Copies of a multi-copy target (rRNA operons) tie and are all kept.
    """
    if not found:
        return []
    fewest = min(p.mismatches_forward + p.mismatches_reverse for p in found)
    return [p for p in found if p.mismatches_forward + p.mismatches_reverse == fewest]


def genome_products(fasta: Path, primers: dict[str, Primer], amplicons: dict[str, Amplicon],
                    max_mismatches: int, max_product: int) -> dict[str, list[Product]]:
    """{amplicon: best products} of one genome. Each primer is searched once per contig, in both orientations."""
    used = {name for a in amplicons.values() for name in (a.forward, a.reverse)}
    found: dict[str, list[Product]] = {name: [] for name in amplicons}
    for contig, seq in read_fasta(fasta).items():
        codes = encode(seq)
        sites: dict[str, Sites] = {name: primer_sites(codes, primers[name].sequence, max_mismatches) for name in sorted(used)}
        for amplicon in amplicons.values():
            found[amplicon.name] += products(contig, seq, sites[amplicon.forward], sites[amplicon.reverse], max_product)
    return {name: best(prods) for name, prods in found.items()}


def product_records(found: dict[str, list[Product]]) -> list[Record]:
    """Distinct best products of one genome as cache records; identical copies are merged and counted."""
    records: list[Record] = []
    for name, prods in found.items():
        merged: dict[str, Record] = {}
        for product in sorted(prods, key=lambda p: (p.contig, p.start)):
            if product.sequence in merged:
                merged[product.sequence]["copies"] = int(merged[product.sequence]["copies"]) + 1
                continue
            merged[product.sequence] = {
                "amplicon": name, "length": product.length, "mismatches_forward": product.mismatches_forward,
                "mismatches_reverse": product.mismatches_reverse, "copies": 1, "contig": product.contig,
                "start": product.start + 1, "end": product.end, "strand": product.strand, "sequence": product.sequence}
        records += list(merged.values())
    return records


def settings_key(primers: dict[str, Primer], amplicons: dict[str, Amplicon], max_mismatches: int, max_product: int) -> str:
    """Changes whenever a primer sequence, a pairing or a search limit changes: cached genomes are then redone."""
    text = json.dumps([[a.name, primers[a.forward].sequence, primers[a.reverse].sequence] for a in amplicons.values()]
                      + [max_mismatches, max_product])
    return hashlib.sha256(text.encode()).hexdigest()[:12]


def load_cache(path: Path, key: str) -> dict[str, list[Record]]:
    """{accession: records} of a per-genome JSONL cache, only entries written with `key`."""
    cache: dict[str, list[Record]] = {}
    if path.is_file():
        with path.open() as handle:
            for line in handle:
                try:
                    entry = json.loads(line)
                except json.JSONDecodeError:  # a line cut off by a killed job
                    continue
                if entry.get("key") == key:
                    cache[entry["accession"]] = entry["records"]
    return cache
