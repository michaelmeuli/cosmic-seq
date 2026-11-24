"""Repo-specific constants; project-wide paths are in `kansasii_paths` (shared, vendored)."""
from pathlib import Path

from .kansasii_paths import GTDB_RELEASE, OUT_COSMIC

REPO = Path(__file__).resolve().parent.parent
PRIMERS_TSV = REPO / "primers.tsv"      # curated: one row per primer
AMPLICONS_TSV = REPO / "amplicons.tsv"  # curated: one row per primer pair
SPECIES_GROUPS_TSV = REPO / "species_groups.tsv"  # curated: groups NCBI lacks, display labels of all groups
EXTRA_GENOMES_TSV = REPO / "extra_genomes.tsv"  # curated: genomes added to the representatives (taxa GTDB merges)
PRIMERS_MD = REPO / "docs" / "PRIMERS.md"  # written by scripts/02

OUTPUT = OUT_COSMIC / "results"
OUTPUT_GTDB = OUTPUT / f"gtdb{GTDB_RELEASE}"  # 01: accession lists
REPRESENTATIVES_TSV = OUTPUT_GTDB / "mycobacteriaceae_representatives.tsv"
OUTPUT_PRIMERS = OUTPUT / "primers"        # 02: check of the primer table against the handbooks
OUTPUT_REFERENCES = OUTPUT / "references"  # 03: reference amplicons
OUTPUT_SPECIES_GROUPS = OUTPUT / "species_groups"  # 04: species -> group (NCBI + curated)
OUTPUT_COSMIC_FILES = OUTPUT / "cosmic_files"  # 05: bundle for cosmic-files res/sequences/references/
H37RV_ACCESSION = "GCF_000195955.2"  # M. tuberculosis H37Rv, the genome the protocols give positions on
