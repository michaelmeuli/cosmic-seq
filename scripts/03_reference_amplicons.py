#!/usr/bin/env python
"""Reference sequences of the Sanger amplicons: in-silico PCR of every primer pair on the GTDB species representatives.

Input: `primers.tsv` / `amplicons.tsv` (checked by 02), the representative list of 01, the curated
`extra_genomes.tsv` (taxa GTDB merges into another species; the label is used as species) and the genomes in the
shared NCBI dataset directory. For each genome and primer pair every product with at most --max-mismatches substitutions
per primer and at most --max-product bp is found; of these the products with the fewest mismatches are kept
(copies of rRNA operons tie; weaker products are paralogs or unrelated loci). Sequences run from the forward
primer to the reverse primer, primers included.

Output directory (default output/cosmic-seq/results/references):
  fasta/<amplicon>.fasta   one record per distinct sequence of a genome: >accession|amplicon|n species [length, mismatches, copies]
  amplicons.tsv            one row per record (amplicon, accession, genus, species, length, mismatches, copies, contig, start, end, strand)
  coverage.tsv             per amplicon: genomes with a product in total, in genus Mycobacterium, without mismatch; length range
  missing_genomes.tsv      representatives without a genome on disk
  per_genome.jsonl         cache: one line per finished genome, so a rerun only processes new genomes

Usage: 03_reference_amplicons.py [--accession ACC ...] [--limit N] [--jobs N] [--out DIR] [--fresh]
Run the full set through `sbatch submit_references.sbatch` (about 20 s per genome and core).
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root: the cosmic_seq package

from cosmic_seq import (AMPLICONS_TSV, EXTRA_GENOMES_TSV, METADATA_TSV, OUTPUT_REFERENCES, PRIMERS_TSV, REPO,  # noqa: E402
                        REPRESENTATIVES_TSV)
from cosmic_seq.genomes import load_genomes, load_metadata  # noqa: E402
from cosmic_seq.insilico import MAX_MISMATCHES, MAX_PRODUCT  # noqa: E402
from cosmic_seq.kansasii_paths import MANIFEST_SCRIPT, NCBI_DATASET  # noqa: E402
from cosmic_seq.tables import (Amplicon, Primer, Record, genome_fasta, genome_products, load_amplicons,  # noqa: E402
                               load_cache, load_primers, product_records, settings_key)

MYCOBACTERIUM = "Mycobacterium"


def run_genome(accession: str, fasta: str, primers: dict[str, Primer], amplicons: dict[str, Amplicon],
               max_mismatches: int, max_product: int) -> tuple[str, list[Record]]:
    """Distinct best products of one genome; identical copies are merged and counted."""
    return accession, product_records(genome_products(Path(fasta), primers, amplicons, max_mismatches, max_product))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--primers", type=Path, default=PRIMERS_TSV)
    parser.add_argument("--amplicons", type=Path, default=AMPLICONS_TSV)
    parser.add_argument("--representatives", type=Path, default=REPRESENTATIVES_TSV, help="accession, genus, species (01)")
    parser.add_argument("--extra", type=Path, default=EXTRA_GENOMES_TSV,
                        help="curated genomes added to the representatives (accession, label, reason); label = species")
    parser.add_argument("--metadata", type=Path, default=METADATA_TSV, help="GTDB metadata rows of step 01")
    parser.add_argument("--genomes", type=Path, default=NCBI_DATASET, help="directory with <accession>/*_genomic.fna")
    parser.add_argument("--accession", nargs="+", help="only these accessions of the representative list")
    parser.add_argument("--limit", type=int, help="only the first N genomes (test runs)")
    parser.add_argument("--max-mismatches", type=int, default=MAX_MISMATCHES, help="substitutions allowed per primer")
    parser.add_argument("--max-product", type=int, default=MAX_PRODUCT, help="longest product (bp)")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--out", type=Path, default=OUTPUT_REFERENCES)
    parser.add_argument("--fresh", action="store_true", help="ignore the cache")
    args = parser.parse_args()

    primers = load_primers(args.primers)
    amplicons = load_amplicons(args.amplicons, primers)
    genomes = [g.as_row() for g in load_genomes(args.representatives, args.extra, load_metadata(args.metadata))]
    if args.accession:
        unknown = sorted(set(args.accession) - {g["accession"] for g in genomes})
        if unknown:
            parser.error(f"not in {args.representatives.name} or {args.extra.name}: {', '.join(unknown)}")
        genomes = [g for g in genomes if g["accession"] in set(args.accession)]
    if args.limit:
        genomes = genomes[:args.limit]

    args.out.mkdir(parents=True, exist_ok=True)
    key = settings_key(primers, amplicons, args.max_mismatches, args.max_product)
    cache_path = args.out / "per_genome.jsonl"
    cache = {} if args.fresh else load_cache(cache_path, key)

    fastas: dict[str, Path] = {}
    missing: list[dict[str, str]] = []
    for genome in genomes:
        fasta = genome_fasta(args.genomes / genome["accession"])
        if fasta is None:
            missing.append(genome)
        else:
            fastas[genome["accession"]] = fasta
    todo = [acc for acc in fastas if acc not in cache]
    print(f"{len(genomes)} genomes, {len(missing)} without a genome on disk, {len(fastas) - len(todo)} cached, {len(todo)} to do", flush=True)

    if todo:
        with cache_path.open("w" if args.fresh else "a") as handle, ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(run_genome, acc, str(fastas[acc]), primers, amplicons, args.max_mismatches, args.max_product)
                       for acc in todo]
            for done, future in enumerate(as_completed(futures), start=1):
                accession, records = future.result()
                cache[accession] = records
                handle.write(json.dumps({"key": key, "accession": accession, "records": records}) + "\n")
                handle.flush()
                if done % 50 == 0 or done == len(todo):
                    print(f"  {done}/{len(todo)} genomes", flush=True)

    # ------------------------------------------------------------------------------------------------ tables
    fasta_dir = args.out / "fasta"
    fasta_dir.mkdir(exist_ok=True)
    for old in fasta_dir.glob("*.fasta"):
        old.unlink()
    columns = ["amplicon", "accession", "genus", "species", "record", "length", "mismatches_forward", "mismatches_reverse",
               "copies", "contig", "start", "end", "strand"]
    by_amplicon: dict[str, list[tuple[dict[str, str], Record]]] = {name: [] for name in amplicons}
    for genome in genomes:
        for record in cache.get(genome["accession"], []):
            if record["amplicon"] in by_amplicon:
                by_amplicon[str(record["amplicon"])].append((genome, record))
    with (args.out / "amplicons.tsv").open("w", newline="") as table, (args.out / "coverage.tsv").open("w", newline="") as cov:
        writer = csv.writer(table, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        coverage = csv.writer(cov, delimiter="\t", lineterminator="\n")
        coverage.writerow(["amplicon", "target", "scope", "documented_bp", "genomes", "genomes_with_product", "mycobacterium_genomes",
                           "mycobacterium_with_product", "with_product_no_mismatch", "distinct_sequences", "min_bp", "median_bp", "max_bp"])
        n_myco = sum(g["genus"] == MYCOBACTERIUM and g["accession"] in fastas for g in genomes)
        for name, amplicon in amplicons.items():
            rows = by_amplicon[name]
            if rows:
                with (fasta_dir / f"{name}.fasta").open("w") as fasta_out:
                    counter: dict[str, int] = {}
                    for genome, record in rows:
                        acc = genome["accession"]
                        counter[acc] = counter.get(acc, 0) + 1
                        record_id = f"{acc}|{name}|{counter[acc]}"
                        fasta_out.write(f">{record_id} {genome['species']} [length={record['length']} "
                                        f"mismatches={record['mismatches_forward']}+{record['mismatches_reverse']} copies={record['copies']}]\n")
                        sequence = str(record["sequence"])
                        fasta_out.write("\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + "\n")
                        writer.writerow([name, acc, genome["genus"], genome["species"], record_id] + [record[c] for c in columns[5:]])
            lengths = sorted(int(record["length"]) for _, record in rows)
            with_product = {g["accession"] for g, _ in rows}
            perfect = {g["accession"] for g, r in rows if r["mismatches_forward"] == 0 and r["mismatches_reverse"] == 0}
            coverage.writerow([name, amplicon.target, amplicon.scope, amplicon.documented_bp or "", len(fastas), len(with_product), n_myco,
                               len({g["accession"] for g, _ in rows if g["genus"] == MYCOBACTERIUM}), len(perfect),
                               len({str(r["sequence"]) for _, r in rows}),
                               lengths[0] if lengths else "", lengths[len(lengths) // 2] if lengths else "", lengths[-1] if lengths else ""])
    with (args.out / "missing_genomes.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["accession", "genus", "species"])
        writer.writerows([g["accession"], g["genus"], g["species"]] for g in missing)

    subprocess.run(["bash", str(MANIFEST_SCRIPT), str(args.out), str(REPO), "scripts/03_reference_amplicons.py",
                    str(args.primers), str(args.amplicons), str(args.representatives), str(args.extra), str(args.genomes)], check=True)
    print(f"{sum(len(rows) for rows in by_amplicon.values())} reference sequences of {len(amplicons)} amplicons in {fasta_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
