"""Genome lists: the GTDB species representatives of step 01, the curated `extra_genomes.tsv`, and the GTDB metadata."""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Genome:
    accession: str
    genus: str
    species: str        # GTDB species, or the label of an extra genome
    extra: bool         # from extra_genomes.tsv (a taxon GTDB merges into another species)
    ncbi_organism: str
    type_strain: bool   # NCBI: assembly from type material

    def as_row(self) -> dict[str, str]:
        return {"accession": self.accession, "genus": self.genus, "species": self.species}


@dataclass(frozen=True)
class Metadata:
    accession: str      # without the GTDB RS_/GB_ prefix
    gtdb_species: str   # "Mycobacterium kansasii"
    ncbi_organism: str
    type_strain: bool
    representative: bool


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_metadata(path: Path) -> dict[str, Metadata]:
    """{accession: metadata} of every genome in the GTDB metadata rows of step 01."""
    metadata: dict[str, Metadata] = {}
    for row in _rows(path):
        accession = row["accession"].split("_", 1)[1] if row["accession"][:3] in ("RS_", "GB_") else row["accession"]
        designation = row.get("ncbi_type_material_designation", "na")
        metadata[accession] = Metadata(accession, row["gtdb_taxonomy"].rsplit(";s__", 1)[-1], row["ncbi_organism_name"],
                                       designation not in ("", "na", "none"), row["gtdb_representative"] == "t")
    return metadata


def load_extra(path: Path) -> list[dict[str, str]]:
    """Rows of `extra_genomes.tsv` (accession, label, reason); empty if the file is absent."""
    if not path.is_file():
        return []
    rows = _rows(path)
    seen: set[str] = set()
    for row in rows:
        if not row["label"].strip() or len(row["label"].split()) < 2:
            raise ValueError(f"{path.name}: {row['accession']}: label must be 'Genus species ...'")
        if row["accession"] in seen:
            raise ValueError(f"{path.name}: {row['accession']} listed twice")
        seen.add(row["accession"])
    return rows


def load_genomes(representatives: Path, extra: Path | None, metadata: dict[str, Metadata]) -> list[Genome]:
    """Representatives of step 01 followed by the extra genomes (their label is the species)."""
    genomes: list[Genome] = []
    for row in _rows(representatives):
        meta = metadata.get(row["accession"])
        genomes.append(Genome(row["accession"], row["genus"], row["species"], False,
                              meta.ncbi_organism if meta else "", meta.type_strain if meta else False))
    known = {g.accession for g in genomes}
    for row in load_extra(extra) if extra else []:
        if row["accession"] in known:
            raise ValueError(f"{row['accession']} of extra_genomes.tsv is already a representative")
        meta = metadata.get(row["accession"])
        label = " ".join(row["label"].split())
        genomes.append(Genome(row["accession"], label.split()[0], label, True,
                              meta.ncbi_organism if meta else "", meta.type_strain if meta else False))
    return genomes


def genome_gff(genome_dir: Path) -> Path | None:
    """The GFF3 annotation of an NCBI datasets directory, if it was delivered."""
    gff = genome_dir / "genomic.gff"
    return gff if gff.is_file() and gff.stat().st_size > 0 else None
