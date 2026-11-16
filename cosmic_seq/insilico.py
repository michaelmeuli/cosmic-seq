"""In-silico PCR on a genome: every product of a primer pair, with the mismatches of both primers.

Started from hsp65-refalign `hsp65_refalign/insilico.py` (itself from mlsa-kansasii `mlsa/loci.py`): same matching
rule (IUPAC codes of the primer honoured, substitutions only, no gaps). Changed here: all products are returned
instead of the first one (16S and 23S have several copies), each primer is searched once per genome in both
orientations and the hits are reused by every amplicon that names the primer, and the search prefilters on the
first positions of the primer so that a family-wide run stays affordable.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt

MAX_MISMATCHES = 3
MAX_PRODUCT = 3000
PREFILTER_LENGTH = 10  # primer positions compared genome-wide before the rest is checked at the surviving sites

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")
_UNKNOWN = 4  # code of every genome base that is not A, C, G or T
_LUT = np.full(256, _UNKNOWN, dtype=np.uint8)
for _code, _base in enumerate("ACGT"):
    _LUT[ord(_base)] = _LUT[ord(_base.lower())] = _code


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def clean_primer(sequence: str) -> str:
    """Primer as written in a protocol (spaces, triplets, lower case) -> upper-case IUPAC string."""
    primer = "".join(sequence.split()).upper().replace("U", "T")
    bad = sorted(set(primer) - set(IUPAC))
    if bad or not primer:
        raise ValueError(f"not a primer sequence: {sequence!r} (unexpected {''.join(bad)!r})")
    return primer


def encode(seq: str) -> npt.NDArray[np.uint8]:
    """Genome sequence as codes 0-3 (A, C, G, T) and 4 (anything else)."""
    return _LUT[np.frombuffer(seq.encode("ascii", errors="replace"), dtype=np.uint8)]


def _mismatch_tables(primer: str) -> npt.NDArray[np.uint8]:
    """[position, genome code] -> 1 if the genome base is not allowed by the primer base at that position."""
    table = np.ones((len(primer), 5), dtype=np.uint8)
    for i, base in enumerate(primer):
        for allowed in IUPAC[base]:
            table[i, "ACGT".index(allowed)] = 0
    return table


def find_sites(codes: npt.NDArray[np.uint8], primer: str, max_mismatches: int = MAX_MISMATCHES) -> dict[int, int]:
    """{0-based start: mismatches} of every window of `codes` that matches `primer` with at most max_mismatches."""
    length = len(primer)
    n = len(codes) - length + 1
    if n <= 0:
        return {}
    table = _mismatch_tables(primer)
    head = min(PREFILTER_LENGTH, length)
    mismatches = np.zeros(n, dtype=np.uint8)
    for i in range(head):
        mismatches += table[i][codes[i:i + n]]
    candidates = np.flatnonzero(mismatches <= max_mismatches)
    counts = mismatches[candidates].astype(np.int64)
    for i in range(head, length):
        counts += table[i][codes[candidates + i]]
    keep = counts <= max_mismatches
    return dict(zip(candidates[keep].tolist(), counts[keep].tolist()))


@dataclass(frozen=True)
class Sites:
    """Binding sites of one primer on one contig, by start on the plus strand."""

    length: int
    plus: dict[int, int]   # primer as written matches the plus strand: it is extended towards higher positions
    minus: dict[int, int]  # reverse complement matches the plus strand: it is extended towards lower positions


def primer_sites(codes: npt.NDArray[np.uint8], primer: str, max_mismatches: int = MAX_MISMATCHES) -> Sites:
    return Sites(len(primer), find_sites(codes, primer, max_mismatches), find_sites(codes, revcomp(primer), max_mismatches))


@dataclass(frozen=True)
class Product:
    contig: str
    start: int   # 0-based start on the plus strand of the contig
    end: int     # exclusive
    strand: str  # "+" if the forward primer binds as written on the plus strand
    mismatches_forward: int
    mismatches_reverse: int
    sequence: str  # 5'->3' from the forward primer to the reverse complement of the reverse primer, primers included

    @property
    def length(self) -> int:
        return self.end - self.start


def products(contig: str, seq: str, forward: Sites, reverse: Sites, max_product: int = MAX_PRODUCT) -> list[Product]:
    """Every product of a primer pair on one contig: each forward site with its nearest reverse site downstream."""
    found: list[Product] = []
    # forward primer on the plus strand, reverse primer on the minus strand
    ends = sorted(start + reverse.length for start in reverse.minus)
    for fstart, fmm in sorted(forward.plus.items()):
        end = next((e for e in ends if e >= fstart + forward.length), None)
        if end is not None and end - fstart <= max_product:
            found.append(Product(contig, fstart, end, "+", fmm, reverse.minus[end - reverse.length], seq[fstart:end].upper()))
    # mirror image: reverse primer on the plus strand, forward primer on the minus strand
    starts = sorted(reverse.plus, reverse=True)
    for fstart, fmm in sorted(forward.minus.items()):
        end = fstart + forward.length
        start = next((s for s in starts if s + reverse.length <= fstart), None)
        if start is not None and end - start <= max_product:
            found.append(Product(contig, start, end, "-", fmm, reverse.plus[start], revcomp(seq[start:end].upper())))
    return found


def read_fasta(path: Path) -> dict[str, str]:
    """{contig id: sequence} of a FASTA file."""
    contigs: dict[str, list[str]] = {}
    current: list[str] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                current = contigs.setdefault(line[1:].split()[0], [])
            else:
                current.append(line.strip())
    return {name: "".join(parts) for name, parts in contigs.items()}
