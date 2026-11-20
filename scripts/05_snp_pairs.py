#!/usr/bin/env python
"""Diagnostic SNP tables: positions that tell two species apart, in the coordinates of the step-04 gene records.

Input: `snp_pairs.tsv` (curated: pair, gene, species A / B, window amplicon, use, species to report, legacy GenBank
records and positions), the step-04 gene records (genes/genes.tsv, genes/fasta), the step-03 products of the window
amplicon, the GTDB metadata of step 01 and the genomes in the shared directory.

Per pair:
  1. The representative gene records of A and B (an extra genome stands for a taxon GTDB merges, e.g. M. ulcerans)
     are aligned; substitutions inside the window (the sequenced amplicon of A; the whole gene for a pair without
     window amplicon) are the candidate positions. A pair with an indel inside the window is reported, not used.
  2. Validation: the gene of **every** genome of A and of B on disk (GTDB species; an extra label matches the NCBI
     organism name, and those genomes are left out of the GTDB species they sit in) is aligned to the A record and the
     alleles at each position are counted. `diagnostic=t` when both species carry their allele in >= --min-fixed of
     the genomes with a call.
  3. Continuity: the legacy GenBank records (fetched once from NCBI into snps/legacy/) are aligned to the A / B
     records and the legacy positions are placed in the new coordinates.

Output directory (default output/cosmic-seq/results/snps):
  pairs.tsv        per pair: records, window, substitutions, indels in the window, diagnostic positions, status
  snp_pairs.tsv    per position: 0-based in A and B, H37Rv gene nt / codon, alleles, diagnostic, alleles of report species
  validation.tsv   per position and species: genomes, calls, counts of A / C / G / T / other
  continuity.tsv   per legacy position: position in A, in the window, substitution, diagnostic
  per_genome.jsonl cache of the gene of non-representative genomes

Usage: 05_snp_pairs.py [--pair NAME ...] [--min-fixed 0.99] [--jobs N] [--out DIR] [--fresh]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
import urllib.request
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root: the cosmic_seq package

from cosmic_seq import (AMPLICONS_TSV, EXTRA_GENOMES_TSV, GENES_TSV, H37RV_ACCESSION, METADATA_TSV,  # noqa: E402
                        OUTPUT_GENES, OUTPUT_REFERENCES, OUTPUT_SNPS, PRIMERS_TSV, REPO, REPRESENTATIVES_TSV,
                        SNP_PAIRS_TSV)
from cosmic_seq.genes import GeneSpec, References, align, genome_genes, load_genes, locus_sequence  # noqa: E402
from cosmic_seq.genomes import Metadata, genome_gff, load_extra, load_metadata  # noqa: E402
from cosmic_seq.insilico import MAX_MISMATCHES, MAX_PRODUCT, read_fasta, revcomp  # noqa: E402
from cosmic_seq.kansasii_paths import MANIFEST_SCRIPT, NCBI_DATASET  # noqa: E402
from cosmic_seq.tables import (Amplicon, Primer, genome_fasta, genome_products, load_amplicons, load_cache,  # noqa: E402
                               load_primers, product_records, settings_key)

VERSION = "1"
Row = list["str | int | float"]
MIN_FIXED = 0.99
BASES = "ACGT"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text&tool=cosmic-seq"


@dataclass(frozen=True)
class Pair:
    name: str
    gene: str
    species_a: str
    species_b: str
    window_amplicon: str  # empty: whole gene
    use: str              # diagnostic / research
    report_species: tuple[str, ...]
    legacy_a: str
    legacy_b: str
    legacy_positions: tuple[int, ...]  # 0-based in legacy_a


def load_pairs(path: Path, genes: set[str], amplicons: set[str]) -> list[Pair]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    pairs: list[Pair] = []
    for row in rows:
        pair = Pair(row["pair"], row["gene"], row["species_a"], row["species_b"], row["window_amplicon"], row["use"],
                    tuple(s.strip() for s in row["report_species"].split(",") if s.strip()), row["legacy_a"],
                    row["legacy_b"], tuple(int(p) for p in row["legacy_positions"].split(",") if p.strip()))
        if pair.gene not in genes:
            raise ValueError(f"{path.name}: {pair.name}: gene {pair.gene!r} not in genes.tsv")
        if pair.window_amplicon and pair.window_amplicon not in amplicons:
            raise ValueError(f"{path.name}: {pair.name}: window amplicon {pair.window_amplicon!r} not in amplicons.tsv")
        if pair.legacy_positions and not pair.legacy_a:
            raise ValueError(f"{path.name}: {pair.name}: legacy positions without legacy_a record")
        pairs.append(pair)
    return pairs


@dataclass(frozen=True)
class GeneRec:
    record_id: str
    accession: str
    sequence: str    # record, 5'->3' in gene orientation
    gene_start: int  # 1-based in the record
    gene_end: int


def load_gene_records(genes_dir: Path) -> dict[str, list[GeneRec]]:
    """{gene: records of step 04}; the order of genes.tsv (record |1 first)."""
    out: dict[str, list[GeneRec]] = {}
    with (genes_dir / "genes.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    fastas: dict[str, dict[str, str]] = {}
    for row in rows:
        if row["gene"] not in fastas:
            fastas[row["gene"]] = read_fasta(genes_dir / "fasta" / f"{row['gene']}.fasta")
        out.setdefault(row["gene"], []).append(GeneRec(row["record"], row["accession"], fastas[row["gene"]][row["record"]],
                                                       int(row["gene_start"]), int(row["gene_end"])))
    return out


def representative_record(species: str, records: list[GeneRec], rows: dict[str, dict[str, str]]) -> GeneRec | None:
    """Record of the species: the GTDB representative, else (for an extra label) the type strain, else the first."""
    candidates = [r for r in records if rows[r.record_id]["species"] == species]
    candidates.sort(key=lambda r: (rows[r.record_id]["extra"] == "t", rows[r.record_id]["type_strain"] != "t",
                                   not r.record_id.endswith("|1")))
    return candidates[0] if candidates else None


def members(species: str, extra_labels: set[str], metadata: dict[str, Metadata], exclude_labels: set[str]) -> list[str]:
    """Accessions of a species: an extra label matches the NCBI organism name, a GTDB species its GTDB taxonomy minus the
    genomes whose NCBI organism is one of `exclude_labels` (the extra taxon of the same pair)."""
    if species in extra_labels:
        return sorted(a for a, m in metadata.items() if m.ncbi_organism == species or m.ncbi_organism.startswith(species + " "))
    return sorted(a for a, m in metadata.items() if m.gtdb_species == species
                  and not any(m.ncbi_organism == x or m.ncbi_organism.startswith(x + " ") for x in exclude_labels))


def genome_gene(accession: str, fasta: str, gff: str | None, spec: GeneSpec, primers: dict[str, Primer],
                amplicons: dict[str, Amplicon], refs: References, max_mismatches: int, max_product: int) -> tuple[str, str]:
    """Gene record of one genome (first record), '' if none: in-silico PCR of the anchors, then step-04 extraction."""
    anchors = {a: amplicons[a] for a in spec.anchors}
    records = product_records(genome_products(Path(fasta), primers, anchors, max_mismatches, max_product))
    genes = [g for g in genome_genes(Path(fasta), Path(gff) if gff else None, {spec.gene: spec}, records, {spec.gene: refs})
             if "failed" not in g]
    return accession, str(genes[0]["sequence"]) if genes else ""


def fetch_legacy(accession: str, directory: Path) -> str:
    path = directory / f"{accession}.fasta"
    if not path.is_file():
        directory.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(EFETCH.format(accession), timeout=60) as response:
            text = response.read().decode()
        if not text.startswith(">"):
            raise SystemExit(f"efetch {accession}: no FASTA returned")
        path.write_text(text)
    return next(iter(read_fasta(path).values())).upper()


def window_in_record(pair: Pair, rec: GeneRec, products: dict[str, list[dict[str, str | int]]]) -> tuple[int, int] | None:
    """0-based [start, end) of the window in the A record."""
    if not pair.window_amplicon:
        return rec.gene_start - 1, rec.gene_end
    found = [p for p in products.get(rec.accession, []) if p["amplicon"] == pair.window_amplicon]
    for product in found:
        for seq in (str(product["sequence"]), revcomp(str(product["sequence"]))):
            pos = rec.sequence.find(seq)
            if pos != -1:
                return pos, pos + len(seq)
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pairs", type=Path, default=SNP_PAIRS_TSV)
    parser.add_argument("--pair", nargs="+", help="only these pairs")
    parser.add_argument("--genes", type=Path, default=GENES_TSV)
    parser.add_argument("--primers", type=Path, default=PRIMERS_TSV)
    parser.add_argument("--amplicons", type=Path, default=AMPLICONS_TSV)
    parser.add_argument("--extra", type=Path, default=EXTRA_GENOMES_TSV)
    parser.add_argument("--metadata", type=Path, default=METADATA_TSV)
    parser.add_argument("--references", type=Path, default=OUTPUT_REFERENCES, help="output directory of step 03")
    parser.add_argument("--gene-records", type=Path, default=OUTPUT_GENES, help="output directory of step 04")
    parser.add_argument("--genomes", type=Path, default=NCBI_DATASET)
    parser.add_argument("--max-mismatches", type=int, default=MAX_MISMATCHES)
    parser.add_argument("--max-product", type=int, default=MAX_PRODUCT)
    parser.add_argument("--min-fixed", type=float, default=MIN_FIXED, help="allele frequency both species must reach")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--out", type=Path, default=OUTPUT_SNPS)
    parser.add_argument("--fresh", action="store_true")
    args = parser.parse_args()

    primers = load_primers(args.primers)
    amplicons = load_amplicons(args.amplicons, primers)
    specs = load_genes(args.genes, set(amplicons))
    pairs = load_pairs(args.pairs, set(specs), set(amplicons))
    if args.pair:
        pairs = [p for p in pairs if p.name in set(args.pair)]
    metadata = load_metadata(args.metadata)
    extra_labels = {" ".join(r["label"].split()) for r in load_extra(args.extra)}
    gene_records = load_gene_records(args.gene_records)
    with (args.gene_records / "genes.tsv").open(newline="") as handle:
        gene_rows = {r["record"]: r for r in csv.DictReader(handle, delimiter="\t")}
    products = load_cache(args.references / "per_genome.jsonl",
                          settings_key(primers, amplicons, args.max_mismatches, args.max_product))
    h37_dir = args.genomes / H37RV_ACCESSION
    h37_fasta, h37_gff = genome_fasta(h37_dir), genome_gff(h37_dir)
    if h37_fasta is None or h37_gff is None:
        raise SystemExit(f"H37Rv genome or GFF missing in {h37_dir}")
    refs = {g: References(locus_sequence(h37_fasta, h37_gff, s.h37rv_locus) if s.h37rv_locus else None,
                          locus_sequence(h37_fasta, h37_gff, s.paralog_locus) if s.paralog_locus else None)
            for g, s in specs.items() if any(p.gene == g for p in pairs)}
    args.out.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------ gene of every genome of the compared species
    step04 = {(r.accession, gene): r.sequence for gene, recs in gene_records.items() for r in reversed(recs)}
    wanted: dict[tuple[str, str], list[str]] = {}  # (pair, species) -> accessions on disk
    for pair in pairs:
        for species, other in ((pair.species_a, pair.species_b), (pair.species_b, pair.species_a)):
            exclude = {other} if other in extra_labels else set()
            wanted[(pair.name, species)] = [a for a in members(species, extra_labels, metadata, exclude)
                                           if genome_fasta(args.genomes / a) is not None]
    cache_path = args.out / "per_genome.jsonl"
    keys = {g: hashlib.sha256(json.dumps([VERSION, asdict(specs[g]), asdict(refs[g]), args.max_mismatches, args.max_product,
                                          [[a, primers[amplicons[a].forward].sequence, primers[amplicons[a].reverse].sequence]
                                           for a in specs[g].anchors]]).encode()).hexdigest()[:12] for g in refs}
    genes_of: dict[tuple[str, str], str] = {}
    if cache_path.is_file() and not args.fresh:
        with cache_path.open() as handle:
            for line in handle:
                try:
                    entry = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if entry["gene"] in keys and entry["key"] == keys[entry["gene"]]:
                    genes_of[(entry["accession"], entry["gene"])] = entry["sequence"]
    todo = sorted({(a, p.gene) for p in pairs for s in (p.species_a, p.species_b) for a in wanted[(p.name, s)]}
                  - set(step04) - set(genes_of))
    print(f"{len(pairs)} pairs; {len(todo)} genomes to extract", flush=True)
    if todo:
        with cache_path.open("w" if args.fresh else "a") as handle, ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = {}
            for acc, gene in todo:
                fasta = genome_fasta(args.genomes / acc)
                gff = genome_gff(args.genomes / acc)
                futures[pool.submit(genome_gene, acc, str(fasta), str(gff) if gff else None, specs[gene], primers,
                                    amplicons, refs[gene], args.max_mismatches, args.max_product)] = gene
            for done, future in enumerate(as_completed(futures), start=1):
                acc, seq = future.result()
                gene = futures[future]
                genes_of[(acc, gene)] = seq
                handle.write(json.dumps({"key": keys[gene], "accession": acc, "gene": gene, "sequence": seq}) + "\n")
                handle.flush()
                if done % 25 == 0 or done == len(todo):
                    print(f"  {done}/{len(todo)} genomes", flush=True)
    genes_of.update(step04)

    # ------------------------------------------------------------------------------------------------ pairs
    pair_rows: list[Row] = []
    snp_rows: list[Row] = []
    val_rows: list[Row] = []
    cont_rows: list[Row] = []
    report_all = sorted({s for p in pairs for s in p.report_species})
    for pair in pairs:
        recs = gene_records.get(pair.gene, [])
        rec_a = representative_record(pair.species_a, recs, gene_rows)
        rec_b = representative_record(pair.species_b, recs, gene_rows)
        n_a, n_b = len(wanted[(pair.name, pair.species_a)]), len(wanted[(pair.name, pair.species_b)])
        if rec_a is None or rec_b is None:
            pair_rows.append([pair.name, pair.gene, pair.use, rec_a.record_id if rec_a else "", rec_b.record_id if rec_b else "",
                              "", "", "", "", 0, n_a, n_b, "no step-04 record for " +
                              " / ".join(s for s, r in ((pair.species_a, rec_a), (pair.species_b, rec_b)) if r is None)])
            continue
        window = window_in_record(pair, rec_a, products)
        if window is None:
            pair_rows.append([pair.name, pair.gene, pair.use, rec_a.record_id, rec_b.record_id, "", "", "", "", 0, n_a, n_b,
                              f"window amplicon {pair.window_amplicon} not found in the A record"])
            continue
        w0, w1 = window
        ab = align(rec_a.sequence, rec_b.sequence).a_to_b
        mapped = [ab[i] for i in range(w0, w1)]
        indels = sum(m is None for m in mapped) + sum(
            1 for x, y in zip(mapped, mapped[1:]) if x is not None and y is not None and y - x != 1)
        subs = [i for i in range(w0, w1) if ab[i] is not None and rec_a.sequence[i] in BASES
                and rec_b.sequence[ab[i] or 0] in BASES and rec_a.sequence[i] != rec_b.sequence[ab[i] or 0]]
        if indels:
            pair_rows.append([pair.name, pair.gene, pair.use, rec_a.record_id, rec_b.record_id, w0, w1, len(subs), indels, 0,
                              n_a, n_b, "indels inside the window: not used"])
            continue
        h37 = refs[pair.gene].gene
        a_h37 = align(rec_a.sequence, h37).a_to_b if h37 else [None] * len(rec_a.sequence)
        report_maps = {}
        for species in pair.report_species:
            rec = representative_record(species, recs, gene_rows)
            if rec is not None:
                report_maps[species] = (rec.sequence, align(rec_a.sequence, rec.sequence).a_to_b)
        # alleles of every genome of A and B at the substitutions
        counts: dict[tuple[str, int], Counter[str]] = {}
        for species in (pair.species_a, pair.species_b):
            for acc in wanted[(pair.name, species)]:
                seq = genes_of.get((acc, pair.gene), "")
                if not seq:
                    for i in subs:
                        counts.setdefault((species, i), Counter())["no_gene"] += 1
                    continue
                to_g = align(rec_a.sequence, seq).a_to_b
                for i in subs:
                    j = to_g[i]
                    base = seq[j] if j is not None else "-"
                    counts.setdefault((species, i), Counter())[base if base in BASES else "other"] += 1
        n_diag = 0
        for i in subs:
            allele_a, allele_b = rec_a.sequence[i], rec_b.sequence[ab[i] or 0]
            fixed = []
            for species, allele in ((pair.species_a, allele_a), (pair.species_b, allele_b)):
                c = counts.get((species, i), Counter())
                called = sum(c[b] for b in BASES)
                fixed.append(c[allele] / called if called else 0.0)
                val_row: Row = [pair.name, i, species, allele, len(wanted[(pair.name, species)]), called]
                val_row += [c[b] for b in BASES]
                val_row += [c["other"] + c["-"], c["no_gene"], round(c[allele] / called, 4) if called else ""]
                val_rows.append(val_row)
            diagnostic = min(fixed) >= args.min_fixed
            n_diag += diagnostic
            h = a_h37[i]
            snp_row: Row = [pair.name, pair.gene, pair.use, rec_a.record_id, rec_b.record_id, i, ab[i] or 0, i - w0,
                             h + 1 if h is not None else "", h // 3 + 1 if h is not None else "",
                             h % 3 + 1 if h is not None else "", allele_a, allele_b, round(fixed[0], 4), round(fixed[1], 4),
                             "t" if diagnostic else "f"]
            snp_row += [(report_maps[s][0][report_maps[s][1][i] or 0] if report_maps[s][1][i] is not None else "-")
                        if s in report_maps else "" for s in report_all]
            snp_rows.append(snp_row)
        pair_rows.append([pair.name, pair.gene, pair.use, rec_a.record_id, rec_b.record_id, w0, w1, len(subs), 0, n_diag,
                          n_a, n_b, "ok"])
        # continuity with the legacy GenBank records
        if pair.legacy_a:
            legacy_a = fetch_legacy(pair.legacy_a, args.out / "legacy")
            to_legacy = align(legacy_a, rec_a.sequence, both_strands=True)  # legacy index -> A index (revcomp if reverse)
            legacy_b = fetch_legacy(pair.legacy_b, args.out / "legacy") if pair.legacy_b else ""
            b_map = align(rec_b.sequence, legacy_b, both_strands=True) if legacy_b else None
            legacy_b_seq = revcomp(legacy_b) if b_map and b_map.reverse else legacy_b
            diag = {int(r[5]): r[15] for r in snp_rows if r[0] == pair.name}
            for pos in pair.legacy_positions:
                j = to_legacy.a_to_b[pos] if pos < len(legacy_a) else None
                if j is not None and to_legacy.reverse:
                    j = len(rec_a.sequence) - 1 - j
                allele_b = ""
                if j is not None and b_map is not None and ab[j] is not None:
                    k = b_map.a_to_b[ab[j] or 0]
                    allele_b = legacy_b_seq[k] if k is not None else "-"
                in_window = j is not None and w0 <= j < w1
                note = ("not aligned" if j is None else "outside the window" if not in_window
                        else "not a substitution between the representatives" if j not in diag
                        else "kept" if diag[j] == "t" else "substitution, not fixed")
                cont_rows.append([pair.name, pair.legacy_a, pos, legacy_a[pos] if pos < len(legacy_a) else "",
                                  pair.legacy_b, allele_b, j if j is not None else "",
                                  rec_a.sequence[j] if j is not None else "", "t" if in_window else "f",
                                  "t" if j in diag else "f", diag.get(j, "") if j is not None else "", note])

    def write(name: str, header: list[str], rows: list[Row]) -> None:
        with (args.out / name).open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            writer.writerows(rows)

    write("pairs.tsv", ["pair", "gene", "use", "record_a", "record_b", "window_start", "window_end", "substitutions",
                        "indels_in_window", "diagnostic", "genomes_a", "genomes_b", "status"], pair_rows)
    write("snp_pairs.tsv", ["pair", "gene", "use", "record_a", "record_b", "pos_a", "pos_b", "pos_window", "h37rv_nt",
                            "h37rv_codon", "h37rv_codon_pos", "allele_a", "allele_b", "fixed_a", "fixed_b", "diagnostic"]
          + [f"allele:{s}" for s in report_all], snp_rows)
    write("validation.tsv", ["pair", "pos_a", "species", "allele", "genomes", "called", "A", "C", "G", "T", "other_or_gap",
                             "no_gene", "fraction"], val_rows)
    write("continuity.tsv", ["pair", "legacy_a", "legacy_pos", "legacy_a_allele", "legacy_b", "legacy_b_allele", "pos_a",
                             "allele_a", "in_window", "substitution", "diagnostic", "note"], cont_rows)
    subprocess.run(["bash", str(MANIFEST_SCRIPT), str(args.out), str(REPO), "scripts/05_snp_pairs.py", str(args.pairs),
                    str(args.gene_records), str(args.references / "per_genome.jsonl"), str(args.metadata), str(args.genomes),
                    "legacy records: NCBI efetch nuccore (snps/legacy/)"], check=True)
    for row in pair_rows:
        print(f"{row[0]}: {row[7]} substitutions in the window, {row[9]} diagnostic ({row[12]})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
