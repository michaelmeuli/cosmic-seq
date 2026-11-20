"""Repo-specific constants; project-wide paths are in `kansasii_paths` (shared, vendored)."""
from pathlib import Path

from .kansasii_paths import GTDB_RELEASE, OUT_COSMIC

REPO = Path(__file__).resolve().parent.parent
PRIMERS_TSV = REPO / "primers.tsv"      # curated: one row per primer
AMPLICONS_TSV = REPO / "amplicons.tsv"  # curated: one row per primer pair
PRIMERS_MD = REPO / "docs" / "PRIMERS.md"  # written by scripts/02
GENES_TSV = REPO / "genes.tsv"              # curated: one row per full-length gene (04)
EXTRA_GENOMES_TSV = REPO / "extra_genomes.tsv"  # curated: genomes added to the representatives (taxa GTDB merges)
SNP_PAIRS_TSV = REPO / "snp_pairs.tsv"      # curated: species pairs told apart by SNPs (05)

OUTPUT = OUT_COSMIC / "results"
OUTPUT_GTDB = OUTPUT / f"gtdb{GTDB_RELEASE}"  # 01: accession lists
REPRESENTATIVES_TSV = OUTPUT_GTDB / "mycobacteriaceae_representatives.tsv"
METADATA_TSV = OUTPUT_GTDB / "mycobacteriaceae_rows_metadata.tsv"  # GTDB metadata of every family genome
OUTPUT_PRIMERS = OUTPUT / "primers"        # 02: check of the primer table against the handbooks
OUTPUT_REFERENCES = OUTPUT / "references"  # 03: reference amplicons
OUTPUT_GENES = OUTPUT / "genes"            # 04: full-length genes
OUTPUT_SNPS = OUTPUT / "snps"              # 05: diagnostic SNP tables
H37RV_ACCESSION = "GCF_000195955.2"  # M. tuberculosis H37Rv, the genome the protocols give positions on
