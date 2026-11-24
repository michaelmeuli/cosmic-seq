"""Genome list of step 03: the GTDB species representatives of step 01 plus the curated `extra_genomes.tsv`."""
from __future__ import annotations

import csv
from pathlib import Path


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_extra(path: Path) -> list[dict[str, str]]:
    """Rows of `extra_genomes.tsv` (accession, label, reason); empty if the file is absent."""
    if not path.is_file():
        return []
    rows = _rows(path)
    seen: set[str] = set()
    for row in rows:
        if len(row["label"].split()) < 2:
            raise ValueError(f"{path.name}: {row['accession']}: label must be 'Genus species ...'")
        if row["accession"] in seen:
            raise ValueError(f"{path.name}: {row['accession']} listed twice")
        seen.add(row["accession"])
    return rows


def load_genomes(representatives: Path, extra: Path | None) -> list[dict[str, str]]:
    """Rows (accession, genus, species) of the representatives, then of the extra genomes (label = species)."""
    genomes = [{"accession": r["accession"], "genus": r["genus"], "species": r["species"]} for r in _rows(representatives)]
    known = {g["accession"] for g in genomes}
    for row in load_extra(extra) if extra else []:
        if row["accession"] in known:
            raise ValueError(f"{row['accession']} of {extra.name if extra else ''} is already a representative")
        label = " ".join(row["label"].split())
        genomes.append({"accession": row["accession"], "genus": label.split()[0], "species": label})
    return genomes
