"""Canonical species names: the genus as GTDB writes it.

GTDB keeps one genus *Mycobacterium*; NCBI uses the 2018 split genera (Gupta et al.), and NCBI deflines also write
`[Mycobacterium] chelonae` for taxa not yet moved. cosmic-files matches species the same way
(`MYCOBACTERIUM_SYNONYMS` in its `src/sequencing/mod.rs`); keep both lists in sync.
"""
from __future__ import annotations

MYCOBACTERIUM_SYNONYMS = ("Mycobacteroides", "Mycolicibacterium", "Mycolicibacter", "Mycolicibacillus")


def canonical_species(name: str) -> str:
    """"Mycobacteroides abscessus" / "[Mycobacterium] chelonae" -> "Mycobacterium abscessus" / "... chelonae"."""
    genus, sep, rest = name.strip().partition(" ")
    genus = genus.strip("[]")
    if genus in MYCOBACTERIUM_SYNONYMS:
        genus = "Mycobacterium"
    return f"{genus}{sep}{rest}"
