#!/usr/bin/env python
"""Full-length gene references: the whole gene around the amplicon products of step 03.

Input: `genes.tsv` (curated: gene, anchor amplicons, feature type, H37Rv locus), the step-03 cache
(references/per_genome.jsonl, same primers and search limits), the genomes and their GFF3 in the shared NCBI dataset
directory, the representatives of 01 and `extra_genomes.tsv`. For each genome and gene the best product of the first
anchor amplicon that has one is taken; the overlapping GFF feature is exported in full (gene orientation), else the
H37Rv gene is aligned around the product (`source=homology`). See cosmic_seq/genes.py.

Output directory (default output/cosmic-seq/results/genes):
  fasta/<gene>.fasta   >accession|gene|n species [gene= length= gene_start= gene_end= source= anchor= copies= type_strain= ncbi_organism=]
  genes.tsv            one row per record (coordinates, offsets, source, identity to H37Rv and to the paralog)
  coverage.tsv         per gene: genomes with an anchor product, with a gene, by source; length range; paralog flags
  missing_genes.tsv    genomes with an anchor product but no gene, and genomes missing from disk / the step-03 cache
  per_genome.jsonl     cache, redone for a genome when genes.tsv, the H37Rv references or its anchor products change

Usage: 04_reference_genes.py [--accession ACC ...] [--genus G] [--limit N] [--jobs N] [--out DIR] [--fresh]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root: the cosmic_seq package

from cosmic_seq import (AMPLICONS_TSV, EXTRA_GENOMES_TSV, GENES_TSV, H37RV_ACCESSION, METADATA_TSV,  # noqa: E402
                        OUTPUT_GENES, OUTPUT_REFERENCES, PRIMERS_TSV, REPO, REPRESENTATIVES_TSV)
from cosmic_seq.genes import GeneRecord, GeneSpec, References, genome_genes, load_genes, locus_sequence  # noqa: E402
from cosmic_seq.genomes import Genome, genome_gff, load_genomes, load_metadata  # noqa: E402
from cosmic_seq.insilico import MAX_MISMATCHES, MAX_PRODUCT  # noqa: E402
from cosmic_seq.kansasii_paths import MANIFEST_SCRIPT, NCBI_DATASET  # noqa: E402
from cosmic_seq.tables import Record, genome_fasta, load_amplicons, load_cache, load_primers, settings_key  # noqa: E402

VERSION = "1"  # bump when the extraction logic changes: every cached genome is redone
MYCOBACTERIUM = "Mycobacterium"
COLUMNS = ["gene", "accession", "record", "genus", "species", "ncbi_organism", "type_strain", "extra", "contig", "start",
           "end", "strand", "length", "gene_start", "gene_end", "gene_bp", "source", "anchor", "copies", "name_ok",
           "identity_h37rv", "identity_paralog", "paralog_flag"]


def h37rv_references(specs: dict[str, GeneSpec], genomes_dir: Path) -> dict[str, References]:
    fasta = genome_fasta(genomes_dir / H37RV_ACCESSION)
    gff = genome_gff(genomes_dir / H37RV_ACCESSION)
    if fasta is None or gff is None:
        raise SystemExit(f"H37Rv {H37RV_ACCESSION}: genome or GFF missing in {genomes_dir}")
    return {g: References(locus_sequence(fasta, gff, s.h37rv_locus) if s.h37rv_locus else None,
                          locus_sequence(fasta, gff, s.paralog_locus) if s.paralog_locus else None)
            for g, s in specs.items()}


def genome_key(specs: dict[str, GeneSpec], refs: dict[str, References], anchors: list[Record], has_gff: bool) -> str:
    text = json.dumps([VERSION, [asdict(s) for s in specs.values()], [asdict(r) for r in refs.values()],
                       sorted(json.dumps(a, sort_keys=True) for a in anchors), has_gff])
    return hashlib.sha256(text.encode()).hexdigest()[:12]


def anchor_records(records: list[Record], specs: dict[str, GeneSpec]) -> list[Record]:
    names = {a for s in specs.values() for a in s.anchors}
    return [r for r in records if r["amplicon"] in names]


def run_genome(accession: str, fasta: str, gff: str | None, specs: dict[str, GeneSpec], anchors: list[Record],
               refs: dict[str, References]) -> tuple[str, list[GeneRecord]]:
    return accession, genome_genes(Path(fasta), Path(gff) if gff else None, specs, anchors, refs)


def header(record_id: str, genome: Genome, r: GeneRecord) -> str:
    gene_bp = int(r["gene_end"]) - int(r["gene_start"]) + 1
    extra = " extra=t" if genome.extra else ""
    return (f">{record_id} {genome.species} [gene={r['gene']} length={r['length']} gene_start={r['gene_start']} "
            f"gene_end={r['gene_end']} gene_bp={gene_bp} source={r['source']} anchor={r['anchor']} copies={r['copies']} "
            f"type_strain={'t' if genome.type_strain else 'f'}{extra} ncbi_organism=\"{genome.ncbi_organism}\"]")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genes", type=Path, default=GENES_TSV)
    parser.add_argument("--primers", type=Path, default=PRIMERS_TSV)
    parser.add_argument("--amplicons", type=Path, default=AMPLICONS_TSV)
    parser.add_argument("--representatives", type=Path, default=REPRESENTATIVES_TSV)
    parser.add_argument("--extra", type=Path, default=EXTRA_GENOMES_TSV)
    parser.add_argument("--metadata", type=Path, default=METADATA_TSV)
    parser.add_argument("--references", type=Path, default=OUTPUT_REFERENCES, help="output directory of step 03")
    parser.add_argument("--genomes", type=Path, default=NCBI_DATASET)
    parser.add_argument("--max-mismatches", type=int, default=MAX_MISMATCHES, help="as used in step 03")
    parser.add_argument("--max-product", type=int, default=MAX_PRODUCT, help="as used in step 03")
    parser.add_argument("--accession", nargs="+")
    parser.add_argument("--genus", help="only this GTDB genus (e.g. Mycobacterium)")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--out", type=Path, default=OUTPUT_GENES)
    parser.add_argument("--fresh", action="store_true", help="ignore the cache")
    args = parser.parse_args()

    primers = load_primers(args.primers)
    amplicons = load_amplicons(args.amplicons, primers)
    specs = load_genes(args.genes, set(amplicons))
    genomes = load_genomes(args.representatives, args.extra, load_metadata(args.metadata))
    if args.accession:
        unknown = sorted(set(args.accession) - {g.accession for g in genomes})
        if unknown:
            parser.error(f"not a representative or extra genome: {', '.join(unknown)}")
        genomes = [g for g in genomes if g.accession in set(args.accession)]
    if args.genus:
        genomes = [g for g in genomes if g.genus == args.genus]
    if args.limit:
        genomes = genomes[:args.limit]

    products = load_cache(args.references / "per_genome.jsonl",
                          settings_key(primers, amplicons, args.max_mismatches, args.max_product))
    refs = h37rv_references(specs, args.genomes)
    args.out.mkdir(parents=True, exist_ok=True)
    cache_path = args.out / "per_genome.jsonl"
    cached: dict[str, tuple[str, list[GeneRecord]]] = {}
    if cache_path.is_file() and not args.fresh:
        with cache_path.open() as handle:
            for line in handle:
                try:
                    entry = json.loads(line)
                except json.JSONDecodeError:
                    continue
                cached[entry["accession"]] = (entry["key"], entry["records"])

    missing: list[tuple[str, Genome, str]] = []
    jobs: dict[str, tuple[Path, Path | None, list[Record], str]] = {}
    results: dict[str, list[GeneRecord]] = {}
    for genome in genomes:
        fasta = genome_fasta(args.genomes / genome.accession)
        if fasta is None:
            missing += [(g, genome, "no genome on disk") for g in specs]
            continue
        if genome.accession not in products:
            missing += [(g, genome, "not in the step-03 cache (run 03)") for g in specs]
            continue
        anchors = anchor_records(products[genome.accession], specs)
        gff = genome_gff(args.genomes / genome.accession)
        key = genome_key(specs, refs, anchors, gff is not None)
        if genome.accession in cached and cached[genome.accession][0] == key:
            results[genome.accession] = cached[genome.accession][1]
        else:
            jobs[genome.accession] = (fasta, gff, anchors, key)
    print(f"{len(genomes)} genomes, {len(results)} cached, {len(jobs)} to do", flush=True)

    if jobs:
        with cache_path.open("w" if args.fresh else "a") as handle, ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(run_genome, acc, str(fasta), str(gff) if gff else None, specs, anchors, refs)
                       for acc, (fasta, gff, anchors, _) in jobs.items()]
            for done, future in enumerate(as_completed(futures), start=1):
                accession, records = future.result()
                results[accession] = records
                handle.write(json.dumps({"key": jobs[accession][3], "accession": accession, "records": records}) + "\n")
                handle.flush()
                if done % 50 == 0 or done == len(jobs):
                    print(f"  {done}/{len(jobs)} genomes", flush=True)

    # ------------------------------------------------------------------------------------------------ tables
    fasta_dir = args.out / "fasta"
    fasta_dir.mkdir(exist_ok=True)
    for old in fasta_dir.glob("*.fasta"):
        old.unlink()
    rows: dict[str, list[tuple[Genome, GeneRecord]]] = {g: [] for g in specs}
    anchored: dict[str, set[str]] = {g: set() for g in specs}
    for genome in genomes:
        for record in results.get(genome.accession, []):
            anchored[str(record["gene"])].add(genome.accession)
            if "failed" in record:
                missing.append((str(record["gene"]), genome, f"anchor {record['failed']} product but no gene"))
            else:
                rows[str(record["gene"])].append((genome, record))
    with (args.out / "genes.tsv").open("w", newline="") as table, (args.out / "coverage.tsv").open("w", newline="") as cov:
        writer = csv.writer(table, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        coverage = csv.writer(cov, delimiter="\t", lineterminator="\n")
        coverage.writerow(["gene", "genomes", "with_anchor_product", "with_gene", "mycobacterium_genomes",
                           "mycobacterium_with_gene", "records", "source_gff", "source_homology", "name_mismatch",
                           "paralog_flag", "anchor_but_no_gene", "distinct_genes", "min_gene_bp", "median_gene_bp",
                           "max_gene_bp"])
        n_on_disk = len(results)
        n_myco = sum(g.genus == MYCOBACTERIUM and g.accession in results for g in genomes)
        for gene in specs:
            with (fasta_dir / f"{gene}.fasta").open("w") as fasta_out:
                counter: dict[str, int] = {}
                for genome, r in rows[gene]:
                    counter[genome.accession] = counter.get(genome.accession, 0) + 1
                    record_id = f"{genome.accession}|{gene}|{counter[genome.accession]}"
                    fasta_out.write(header(record_id, genome, r) + "\n")
                    seq = str(r["sequence"])
                    fasta_out.write("\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n")
                    gene_bp = int(r["gene_end"]) - int(r["gene_start"]) + 1
                    writer.writerow([gene, genome.accession, record_id, genome.genus, genome.species, genome.ncbi_organism,
                                     "t" if genome.type_strain else "f", "t" if genome.extra else "f", r["contig"],
                                     r["start"], r["end"], r["strand"], r["length"], r["gene_start"], r["gene_end"],
                                     gene_bp, r["source"], r["anchor"], r["copies"], r["name_ok"], r["identity_h37rv"],
                                     r["identity_paralog"], r["paralog_flag"]])
            with_gene = {g.accession for g, _ in rows[gene]}
            bps = sorted(int(r["gene_end"]) - int(r["gene_start"]) + 1 for _, r in rows[gene])
            genes_seq = {str(r["sequence"])[int(r["gene_start"]) - 1:int(r["gene_end"])] for _, r in rows[gene]}
            coverage.writerow([gene, n_on_disk, len(anchored[gene]), len(with_gene), n_myco,
                               len({g.accession for g, _ in rows[gene] if g.genus == MYCOBACTERIUM}), len(rows[gene]),
                               sum(r["source"] == "gff" for _, r in rows[gene]),
                               sum(r["source"] == "homology" for _, r in rows[gene]),
                               sum(r["name_ok"] == "f" for _, r in rows[gene]),
                               sum(r["paralog_flag"] == "t" for _, r in rows[gene]),
                               len(anchored[gene] - with_gene), len(genes_seq),
                               bps[0] if bps else "", bps[len(bps) // 2] if bps else "", bps[-1] if bps else ""])
    with (args.out / "missing_genes.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene", "accession", "genus", "species", "reason"])
        writer.writerows([gene, g.accession, g.genus, g.species, reason] for gene, g, reason in missing)

    subprocess.run(["bash", str(MANIFEST_SCRIPT), str(args.out), str(REPO), "scripts/04_reference_genes.py",
                    str(args.genes), str(args.references / "per_genome.jsonl"), str(args.representatives),
                    str(args.extra), str(args.genomes)], check=True)
    print(f"{sum(len(r) for r in rows.values())} gene records of {len(specs)} genes in {fasta_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
