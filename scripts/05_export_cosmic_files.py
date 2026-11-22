#!/usr/bin/env python
"""Bundle for cosmic-files: the files of its `res/sequences/references/`, in one flat directory.

- `<amplicon>.fasta`: the reference amplicons of step 3 (`references/fasta/`);
- `species_groups.tsv`: the species -> group lookup of step 4 (cosmic-files names a 16S call by its group);
- `MANIFEST.txt`.

FASTA files that step 3 no longer writes are removed from the bundle, so the bundle mirrors the results. Needs no genome
data: runs on the cluster after 03/04, or on Windows with `--results` pointing at a copy of `output/cosmic-seq/results`.
cosmic-files copies the bundle with its `scripts/sync_cosmic_seq.py`.

Usage: 05_export_cosmic_files.py [--results DIR] [--out DIR]
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root: the cosmic_seq package

from cosmic_seq import OUTPUT, OUTPUT_COSMIC_FILES, OUTPUT_REFERENCES, OUTPUT_SPECIES_GROUPS  # noqa: E402
from cosmic_seq.manifest import write_manifest  # noqa: E402

SPECIES_GROUPS_HEADER = ["group", "label", "species", "source", "group_taxid", "species_taxid", "ncbi_assemblies",
                         "ncbi_name"]


def check_fasta(path: Path) -> int:
    """Number of records; fails on an empty file or a header that is not `ACC|amplicon|N Genus species ...`."""
    records = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            first, _, rest = line[1:].partition(" ")
            if first.count("|") != 2 or len(rest.split()) < 2:
                raise ValueError(f"{path}: unexpected header {line.strip()!r}")
            records += 1
    if records == 0:
        raise ValueError(f"{path}: no records")
    return records


def check_species_groups(path: Path) -> int:
    """Number of species; fails on another header or a species listed twice."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].split("\t") != SPECIES_GROUPS_HEADER:
        raise ValueError(f"{path}: header is not {SPECIES_GROUPS_HEADER} (rerun scripts/04)")
    species = [line.split("\t")[2] for line in lines[1:] if line]
    duplicates = sorted({name for name in species if species.count(name) > 1})
    if duplicates:
        raise ValueError(f"{path}: species in several groups: {', '.join(duplicates)}")
    return len(species)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", type=Path, default=OUTPUT,
                        help="results directory of cosmic-seq (with references/ and species_groups/)")
    parser.add_argument("--out", type=Path, help=f"bundle directory (default: {OUTPUT_COSMIC_FILES}, "
                        "or RESULTS/cosmic_files when --results is given)")
    args = parser.parse_args()
    fasta_dir = args.results / OUTPUT_REFERENCES.name / "fasta"
    groups_tsv = args.results / OUTPUT_SPECIES_GROUPS.name / "species_groups.tsv"
    out: Path = args.out or (args.results / OUTPUT_COSMIC_FILES.name)

    sources = sorted(fasta_dir.glob("*.fasta"))
    if not sources:
        raise FileNotFoundError(f"no FASTA in {fasta_dir} (run scripts/03)")
    records = {path.name: check_fasta(path) for path in sources}
    n_species = check_species_groups(groups_tsv)

    out.mkdir(parents=True, exist_ok=True)
    stale = sorted(path for path in out.glob("*.fasta") if path.name not in records)
    for path in stale:
        path.unlink()
        print(f"removed {path.name}")
    for path in sources:
        shutil.copyfile(path, out / path.name)
    shutil.copyfile(groups_tsv, out / groups_tsv.name)

    write_manifest(out, "scripts/05_export_cosmic_files.py", str(fasta_dir), str(groups_tsv))
    print(f"{len(sources)} FASTA ({sum(records.values())} records), species_groups.tsv ({n_species} species) in {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
