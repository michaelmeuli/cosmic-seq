#!/usr/bin/env python
"""Species groups (complexes) of the family Mycobacteriaceae: NCBI taxonomy plus the curated table of the repo.

1. NCBI Datasets taxonomy API: every SPECIES_GROUP / SPECIES_SUBGROUP node below the family (found, not listed by hand);
2. for each group, every node of rank SPECIES below it; placeholders ("Mycobacterium sp. ...", "cf.", quoted isolate
   names) are counted in groups.tsv but left out of species_groups.tsv. Each SUBSPECIES "G e subsp. s" below it is
   also listed as species "G s" unless that species exists: NCBI demoted former species to subspecies (M. chimaera,
   M. yongonense -> M. intracellulare subsp.; M. massiliense, M. bolletii -> M. abscessus subsp.), and old deflines
   still use the species names;
3. NCBI taxa of another rank that serve as a complex, named by a `taxid` row of the curated table (genus
   Mycobacteroides, 670516 = the M. chelonae-abscessus complex): their species are collected the same way (source
   ncbi); the node's NCBI name must equal the row's `group`;
4. the curated `species_groups.tsv` of the repo adds the groups NCBI lacks (kansasii) and synonyms, and gives every
   group its display label. A curated species that is also in an NCBI group is an error unless its note says
   `override`.

Species are written with the genus as GTDB has it (`cosmic_seq.names.canonical_species`: Mycobacteroides,
Mycolicibacterium, ... -> Mycobacterium); column `ncbi_name` keeps the NCBI spelling where it differs.

One row per species and one group per species: the table is a lookup "species name -> group" (cosmic-files names a
16S call by its group, as 16S does not separate the members of a complex). GTDB is not used: it has no rank for
complexes, and the NCBI groups are NCBI names.

Output directory: species_groups.tsv (one row per species), groups.tsv (one row per group), raw/<taxid>.json (the API
answers; --from-raw parses them again without network), MANIFEST.txt (skipped where bash or the manifest script is
missing, e.g. on Windows).

Usage: 04_ncbi_species_groups.py [--family-taxid N] [--curated FILE] [--out DIR] [--from-raw]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root: the cosmic_seq package

from cosmic_seq import OUTPUT_SPECIES_GROUPS, SPECIES_GROUPS_TSV  # noqa: E402
from cosmic_seq.manifest import write_manifest  # noqa: E402
from cosmic_seq.names import canonical_species  # noqa: E402

MYCOBACTERIACEAE_TAXID = 1762
API = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon"
USER_AGENT = "cosmic-seq/1.0 (python urllib)"
GROUP_RANKS = ("SPECIES_GROUP", "SPECIES_SUBGROUP")
RETRIES = 3

Edges = dict[str, dict[str, Any]]


@dataclass
class Member:
    group: str
    label: str
    species: str
    source: str  # ncbi | curated
    group_taxid: str = ""
    species_taxid: str = ""
    ncbi_assemblies: str = ""
    ncbi_name: str = ""  # NCBI spelling when it differs from `species`


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
        data: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    else:
        query = "&".join(f"rank_limits={rank}" for rank in ranks)
        data = fetch(f"{API}/{taxid}/filtered_subtree?{query}")
        part = path.with_suffix(".json.part")
        part.write_text(json.dumps(data, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        part.replace(path)
    edges: Edges = data.get("edges", {})
    if not edges:
        raise RuntimeError(f"taxon {taxid}: no subtree in the answer ({path})")
    return edges


def is_formal(name: str) -> bool:
    """A validly named species, not a placeholder ("Mycobacterium sp. X", "cf. triplex", "'isolate 23'")."""
    return " sp." not in name and " cf." not in name and "'" not in name and '"' not in name


def subspecies_as_species(name: str) -> str | None:
    """"Mycobacterium intracellulare subsp. chimaera" -> "Mycobacterium chimaera" (canonical genus); None if no subsp."""
    head, sep, rest = name.partition(" subsp. ")
    words = rest.split()
    if not sep or not words:
        return None
    genus = canonical_species(head).split()[0]
    return f"{genus} {words[0]}"


def short_label(group: str) -> str:
    """Display label of a group without a curated one: "Mycobacterium simiae complex" -> "M. simiae complex"."""
    genus, _, rest = group.partition(" ")
    return f"{genus[0]}. {rest}" if rest else group


def read_curated(path: Path) -> tuple[dict[str, str], list[tuple[Member, str]], dict[int, str]]:
    """Curated table: group -> label, the curated species with their note, and taxid -> group of extra NCBI taxa."""
    labels: dict[str, str] = {}
    members: list[tuple[Member, str]] = []
    extra_taxa: dict[int, str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for line, row in enumerate(csv.DictReader(handle, delimiter="\t"), start=2):
            group, label, species = row["group"].strip(), row["label"].strip(), row["species"].strip()
            taxid = (row.get("taxid") or "").strip()
            if not group or not label:
                raise ValueError(f"{path}:{line}: group and label are required")
            if labels.setdefault(group, label) != label:
                raise ValueError(f"{path}:{line}: group {group!r} has labels {labels[group]!r} and {label!r}")
            if taxid and species:
                raise ValueError(f"{path}:{line}: a row has either a species or a taxid")
            if taxid:
                if int(taxid) in extra_taxa or group in extra_taxa.values():
                    raise ValueError(f"{path}:{line}: taxid {taxid} / group {group!r} listed twice")
                extra_taxa[int(taxid)] = group
            if species:
                members.append((Member(group, label, canonical_species(species), "curated"), row.get("note") or ""))
    return labels, members, extra_taxa


def collect_group(group_taxid: int, group_name: str, label: str, edges: Edges,
                  by_species: dict[str, Member], subspecies_of: Callable[[int], Edges]) -> list[object]:
    """Add the formal species below an NCBI taxon to `by_species`, and their subspecies under the species name the
    subspecies epithet makes (`subspecies_as_species`); returns the groups.tsv row. `subspecies_of(taxid)` answers
    the subtree of a species: the group subtree stops at rank SPECIES."""
    species = [(taxid, edge) for taxid, edge in edges.items() if edge.get("rank") == "SPECIES"]
    formal = [(taxid, edge) for taxid, edge in species if is_formal(str(edge["scientific_name"]))]
    for taxid, edge in formal:
        ncbi_name = str(edge["scientific_name"])
        name = canonical_species(ncbi_name)
        if name in by_species:
            raise ValueError(f"{name}: in NCBI groups {by_species[name].group!r} and {group_name!r}")
        by_species[name] = Member(group_name, label, name, "ncbi", str(group_taxid), taxid,
                                  str(edge.get("assembly_count", 0)), ncbi_name if ncbi_name != name else "")
    subspecies: list[tuple[str, Member]] = []
    for taxid, edge in formal:
        if edge.get("children_status") == "NO_VISIBLE_CHILDREN":
            continue
        for sub_taxid, sub_edge in sorted(subspecies_of(int(taxid)).items()):
            ncbi_name = str(sub_edge["scientific_name"])
            name = subspecies_as_species(ncbi_name) if sub_edge.get("rank") == "SUBSPECIES" else None
            if name is not None and is_formal(ncbi_name) and name != canonical_species(str(edge["scientific_name"])):
                subspecies.append((name, Member(group_name, label, name, "ncbi", str(group_taxid), sub_taxid,
                                                str(sub_edge.get("assembly_count", 0)), ncbi_name)))
    n_subspecies = 0
    for name, member in subspecies:
        if name in by_species:
            continue  # a species of that name (or an earlier subspecies) wins
        by_species[name] = member
        n_subspecies += 1
    n_assembly = sum(bool(edge.get("assembly_count")) for _, edge in formal)
    print(f"{group_name} ({group_taxid}): {len(species)} species, {len(formal)} formal "
          f"({len(species) - len(formal)} placeholders left out), {n_assembly} with assemblies, "
          f"{n_subspecies} subspecies listed as species")
    return [group_name, label, "ncbi", group_taxid, len(species), len(formal), n_assembly]


def write_tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    part = path.with_suffix(path.suffix + ".part")
    with part.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    part.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--family-taxid", type=int, default=MYCOBACTERIACEAE_TAXID)
    parser.add_argument("--curated", type=Path, default=SPECIES_GROUPS_TSV)
    parser.add_argument("--out", type=Path, default=OUTPUT_SPECIES_GROUPS)
    parser.add_argument("--from-raw", action="store_true", help="parse the saved API answers in OUT/raw, no network")
    args = parser.parse_args()

    raw_dir = args.out / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    species_raw_dir = raw_dir / "species"  # subtrees of species with children (their subspecies)
    species_raw_dir.mkdir(exist_ok=True)

    def subspecies_of(taxid: int) -> Edges:
        return subtree(taxid, ("SUBSPECIES",), species_raw_dir, args.from_raw)

    labels, curated, extra_taxa = read_curated(args.curated)

    family = subtree(args.family_taxid, GROUP_RANKS, raw_dir, args.from_raw)
    groups = sorted(((int(taxid), edge) for taxid, edge in family.items() if edge.get("rank") in GROUP_RANKS),
                    key=lambda item: str(item[1]["scientific_name"]))
    if not groups:
        raise RuntimeError(f"no {'/'.join(GROUP_RANKS)} node below taxon {args.family_taxid}")

    by_species: dict[str, Member] = {}
    group_rows: list[list[object]] = []
    for group_taxid, group in groups:
        group_name = str(group["scientific_name"])
        label = labels.get(group_name) or short_label(group_name)
        edges = subtree(group_taxid, ("SPECIES",), raw_dir, args.from_raw)
        group_rows.append(collect_group(group_taxid, group_name, label, edges, by_species, subspecies_of))

    ncbi_groups = {str(group["scientific_name"]) for _, group in groups}
    for taxid, group_name in sorted(extra_taxa.items()):
        if taxid in {group_taxid for group_taxid, _ in groups}:
            raise ValueError(f"taxid {taxid} ({group_name!r}) is already an NCBI species group; drop its taxid")
        edges = subtree(taxid, ("SPECIES",), raw_dir, args.from_raw)
        node_name = str(edges.get(str(taxid), {}).get("scientific_name", ""))
        if node_name != group_name:
            raise ValueError(f"taxid {taxid} is {node_name!r} in NCBI, not {group_name!r} as in {args.curated}")
        group_rows.append(collect_group(taxid, group_name, labels[group_name], edges, by_species,
                                        subspecies_of))
        ncbi_groups.add(group_name)
    for group_name in sorted(set(labels) - ncbi_groups - {member.group for member, _ in curated}):
        print(f"warning: curated label for {group_name!r}, which is no NCBI group and has no curated species",
              file=sys.stderr)

    curated_counts: dict[str, int] = {}
    for member, note in curated:
        known = by_species.get(member.species)
        if known is not None and known.source == "curated":
            raise ValueError(f"{member.species}: listed twice in the curated table (genus synonyms count as one)")
        if known is not None and "override" not in note:
            raise ValueError(f"{member.species}: curated in {member.group!r} but already in {known.group!r} "
                             f"({known.source}); write `override` in the note to move it")
        if known is not None:
            member.ncbi_name = known.ncbi_name
        by_species[member.species] = member
        if member.group not in ncbi_groups:
            curated_counts[member.group] = curated_counts.get(member.group, 0) + 1
    for group_name, count in sorted(curated_counts.items()):
        group_rows.append([group_name, labels[group_name], "curated", "", count, count, ""])
        print(f"{group_name} (curated): {count} species")

    rows = sorted(by_species.values(), key=lambda m: (m.group, m.species))
    write_tsv(args.out / "species_groups.tsv",
              ["group", "label", "species", "source", "group_taxid", "species_taxid", "ncbi_assemblies", "ncbi_name"],
              [[m.group, m.label, m.species, m.source, m.group_taxid, m.species_taxid, m.ncbi_assemblies, m.ncbi_name]
               for m in rows])
    write_tsv(args.out / "groups.tsv",
              ["group", "label", "source", "group_taxid", "species", "formal_species", "species_with_assemblies"],
              group_rows)

    write_manifest(args.out, "scripts/04_ncbi_species_groups.py", str(args.curated),
                   f"{API}/{args.family_taxid}/filtered_subtree")
    print(f"{len(rows)} species of {len(group_rows)} groups in {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
