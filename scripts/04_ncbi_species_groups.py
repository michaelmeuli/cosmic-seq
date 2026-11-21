#!/usr/bin/env python
"""Species lists of the NCBI species groups (complexes) of the family Mycobacteriaceae, mapped to GTDB.

1. NCBI Datasets taxonomy API: every SPECIES_GROUP / SPECIES_SUBGROUP node below the family (found, not listed by hand);
2. for each group, every node of rank SPECIES below it (placeholders such as "Mycobacterium sp. ..." included,
   marked by the column `formal`);
3. each NCBI species is looked up in the GTDB metadata by `ncbi_species_taxid`: how many genomes of the family carry
   it, in which GTDB species they are, and the GTDB species representatives of those species.

NCBI groups are NCBI taxonomy: one GTDB species can hold several NCBI species (all of the M. tuberculosis complex is
GTDB "Mycobacterium tuberculosis", the M. ulcerans group is GTDB "Mycobacterium marinum"). Genomes filed under a wrong
NCBI taxon appear as odd GTDB species with small counts; they are shown, not filtered.

Output directory: species_groups.tsv (one row per group and species), groups.tsv (one row per group),
raw/<taxid>.json (the API answers; --from-raw parses them again without network), MANIFEST.txt.

Usage: 04_ncbi_species_groups.py [--family-taxid N] [--metadata FILE] [--out DIR] [--from-raw]
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
import urllib.error
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root: the cosmic_seq package

from cosmic_seq import OUTPUT_SPECIES_GROUPS, REPO  # noqa: E402
from cosmic_seq.kansasii_paths import GTDB_METADATA_DIR, GTDB_RELEASE, MANIFEST_SCRIPT  # noqa: E402

MYCOBACTERIACEAE_TAXID = 1762
API = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon"
USER_AGENT = "cosmic-seq/1.0 (python urllib)"
GROUP_RANKS = ("SPECIES_GROUP", "SPECIES_SUBGROUP")
RETRIES = 3
FAMILY = "f__Mycobacteriaceae"

Edges = dict[str, dict[str, Any]]


def fetch(url: str) -> dict[str, Any]:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    for attempt in range(1, RETRIES + 1):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                data: dict[str, Any] = json.load(response)
                return data
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as error:
            if attempt == RETRIES:
                raise RuntimeError(f"{url}: {error}") from error
            print(f"  {url}: {error}; retry {attempt}/{RETRIES - 1}", file=sys.stderr)
            time.sleep(5 * attempt)
    raise AssertionError("unreachable")


def subtree(taxid: int, ranks: tuple[str, ...], raw_dir: Path, from_raw: bool) -> Edges:
    """`edges` of the filtered subtree of a taxon, from the API (answer saved to raw_dir) or from raw_dir."""
    path = raw_dir / f"{taxid}.json"
    if from_raw:
        data: dict[str, Any] = json.loads(path.read_text())
    else:
        query = "&".join(f"rank_limits={rank}" for rank in ranks)
        data = fetch(f"{API}/{taxid}/filtered_subtree?{query}")
        part = path.with_suffix(".json.part")
        part.write_text(json.dumps(data, indent=1, sort_keys=True) + "\n")
        part.replace(path)
    edges: Edges = data.get("edges", {})
    if not edges:
        raise RuntimeError(f"taxon {taxid}: no subtree in the answer ({path})")
    return edges


def is_formal(name: str) -> bool:
    """A validly named species, not a placeholder ("Mycobacterium sp. X", "cf. triplex", "'isolate 23'")."""
    return " sp." not in name and " cf." not in name and "'" not in name and '"' not in name


def gtdb_by_species_taxid(metadata: Path) -> tuple[dict[str, Counter[str]], dict[str, str]]:
    """ncbi_species_taxid -> Counter of GTDB species of the family; GTDB species -> its representative accession."""
    by_taxid: dict[str, Counter[str]] = {}
    representative: dict[str, str] = {}
    with metadata.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            taxonomy = row["gtdb_taxonomy"]
            if FAMILY not in taxonomy:
                continue
            species = taxonomy.rsplit(";", 1)[-1].removeprefix("s__")
            by_taxid.setdefault(row["ncbi_species_taxid"], Counter())[species] += 1
            if row["gtdb_representative"] == "t":
                representative[species] = row["accession"][3:]  # without the RS_/GB_ flag of GTDB
    return by_taxid, representative


def write_tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    part = path.with_suffix(path.suffix + ".part")
    with part.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    part.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--family-taxid", type=int, default=MYCOBACTERIACEAE_TAXID)
    parser.add_argument("--metadata", type=Path, default=GTDB_METADATA_DIR / f"bac120_metadata_r{GTDB_RELEASE}.tsv")
    parser.add_argument("--out", type=Path, default=OUTPUT_SPECIES_GROUPS)
    parser.add_argument("--from-raw", action="store_true", help="parse the saved API answers in OUT/raw, no network")
    args = parser.parse_args()

    raw_dir = args.out / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    family = subtree(args.family_taxid, GROUP_RANKS, raw_dir, args.from_raw)
    groups = sorted(((int(taxid), edge) for taxid, edge in family.items() if edge.get("rank") in GROUP_RANKS),
                    key=lambda item: str(item[1]["scientific_name"]))
    if not groups:
        raise RuntimeError(f"no {'/'.join(GROUP_RANKS)} node below taxon {args.family_taxid}")

    by_taxid, representative = gtdb_by_species_taxid(args.metadata)

    species_rows: list[list[object]] = []
    group_rows: list[list[object]] = []
    for group_taxid, group in groups:
        edges = subtree(group_taxid, ("SPECIES",), raw_dir, args.from_raw)
        species = sorted(((taxid, edge) for taxid, edge in edges.items() if edge.get("rank") == "SPECIES"),
                         key=lambda item: str(item[1]["scientific_name"]))
        gtdb_hit: set[str] = set()
        for taxid, edge in species:
            name = str(edge["scientific_name"])
            gtdb = by_taxid.get(taxid, Counter())
            ranked = sorted(gtdb.items(), key=lambda item: (-item[1], item[0]))
            gtdb_hit.update(gtdb)
            species_rows.append([group_taxid, group["scientific_name"], group["rank"], taxid, name,
                                 "yes" if is_formal(name) else "no", edge.get("assembly_count", 0), sum(gtdb.values()),
                                 "; ".join(f"{sp} ({n})" for sp, n in ranked),
                                 "; ".join(representative.get(sp, "") for sp, _ in ranked)])
        n_formal = sum(is_formal(str(edge["scientific_name"])) for _, edge in species)
        n_assembly = sum(bool(edge.get("assembly_count")) for _, edge in species)
        group_rows.append([group_taxid, group["scientific_name"], group["rank"], len(species), n_formal, n_assembly,
                           len(gtdb_hit), "; ".join(sorted(gtdb_hit))])
        print(f"{group['scientific_name']} ({group_taxid}): {len(species)} species, {n_formal} formal, "
              f"{n_assembly} with assemblies, {len(gtdb_hit)} GTDB species")

    write_tsv(args.out / "species_groups.tsv",
              ["group_taxid", "group_name", "group_rank", "species_taxid", "species_name", "formal", "ncbi_assemblies",
               "gtdb_genomes", "gtdb_species", "gtdb_representatives"], species_rows)
    write_tsv(args.out / "groups.tsv",
              ["group_taxid", "group_name", "group_rank", "species", "formal_species", "species_with_assemblies",
               "gtdb_species_count", "gtdb_species"], group_rows)

    subprocess.run(["bash", str(MANIFEST_SCRIPT), str(args.out), str(REPO), "scripts/04_ncbi_species_groups.py",
                    str(args.metadata), f"{API}/{args.family_taxid}/filtered_subtree"], check=True)
    print(f"{len(species_rows)} species rows of {len(groups)} groups in {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
