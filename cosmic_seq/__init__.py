"""Repo-specific constants; project-wide paths are in `kansasii_paths` (shared, vendored)."""
from pathlib import Path

from .kansasii_paths import GTDB_RELEASE, OUT_COSMIC

REPO = Path(__file__).resolve().parent.parent
PRIMERS_TSV = REPO / "primers.tsv"      # curated: one row per primer
AMPLICONS_TSV = REPO / "amplicons.tsv"  # curated: one row per primer pair
PRIMERS_MD = REPO / "docs" / "PRIMERS.md"  # written by scripts/02

OUTPUT = OUT_COSMIC / "results"
OUTPUT_GTDB = OUTPUT / f"gtdb{GTDB_RELEASE}"  # 01: accession lists
REPRESENTATIVES_TSV = OUTPUT_GTDB / "mycobacteriaceae_representatives.tsv"
OUTPUT_PRIMERS = OUTPUT / "primers"        # 02: check of the primer table against the handbooks
OUTPUT_REFERENCES = OUTPUT / "references"  # 03: reference amplicons
OUTPUT_SPECIES_GROUPS = OUTPUT / "species_groups"  # 04: species of the NCBI species groups
H37RV_ACCESSION = "GCF_000195955.2"  # M. tuberculosis H37Rv, the genome the protocols give positions on
