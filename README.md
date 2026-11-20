# cosmic-seq

Reference sequences for the mycobacterial Sanger assays of the lab: the primers documented in the lab handbooks
(`data/Arbeitshandbücher`) are run in silico on the GTDB species representatives of the family Mycobacteriaceae, and
each product is kept as a reference for the reads of that assay.

Five steps:

| Step | Script | Does | Output |
|---|---|---|---|
| 1 | `scripts/01_gtdb_mycobacteriaceae_genomes.sh` | GTDB r232 metadata -> all `f__Mycobacteriaceae` accessions -> NCBI `datasets`, **only genomes that are not on disk** | genomes in the shared `data/gtdb_genomes/Mycobacteriaceae/ncbi_dataset/data/`; lists in `output/cosmic-seq/results/gtdb232/` |
| 2 | `scripts/02_handbook_primers.py` | checks `primers.tsv` against the protocol files, lists other oligos of the mycobacterial protocols, runs every pair on H37Rv | `docs/PRIMERS.md`; `output/cosmic-seq/results/primers/` |
| 3 | `scripts/03_reference_amplicons.py` | in-silico PCR of every primer pair on every species representative (+ `extra_genomes.tsv`) | `output/cosmic-seq/results/references/` (`fasta/<amplicon>.fasta`, `amplicons.tsv`, `coverage.tsv`) |
| 4 | `scripts/04_reference_genes.py` | the full-length gene around the products of its anchor amplicons (`genes.tsv`) | `output/cosmic-seq/results/genes/` (`fasta/<gene>.fasta`, `genes.tsv`, `coverage.tsv`, `missing_genes.tsv`) |
| 5 | `scripts/05_snp_pairs.py` | SNPs that tell two species apart (`snp_pairs.tsv`), validated on every genome of both species | `output/cosmic-seq/results/snps/` (`pairs.tsv`, `snp_pairs.tsv`, `validation.tsv`, `continuity.tsv`) |

## Primers

[`docs/PRIMERS.md`](docs/PRIMERS.md) is the list (primers, primer pairs, protocol of each, product on H37Rv). It is written
by step 2 from two curated tables:

- `primers.tsv`: one row per primer (name, target, scope, direction, sequence, role, protocol, protocol file);
- `amplicons.tsv`: one row per primer pair (forward, reverse, product length of the protocol, use, sequencing primers).

The tables were written by reading the protocols; the notation differs too much between them to parse blindly. Step 2
fails if a sequence of `primers.tsv` is not in the text of its protocol file, so the tables cannot drift from the
handbook unnoticed. Scope: assays on mycobacteria that end in Sanger sequencing (chapters 32.433xx, 32.435xx, 32.436xx of
the molecular handbook, and the research protocols kept outside it, labelled `NZM <n>`). Probes, real-time PCR and
fragment-length assays are listed in `docs/PRIMERS.md` but give no reference.

To add an assay: add its rows to both tables, run step 2, then step 3.

## Genomes

Step 1 follows `immensekansasii/scripts/01_gtdb_mycobacteria_genomes.sh` and fills the same directory, so the 18,613
genomes downloaded for the kansasii project are reused. A genome counts as present when its directory holds a non-empty
`*_genomic.fna`. Accessions that NCBI no longer delivers are written to `unavailable_accessions.txt` and not requested
again (`--retry-unavailable` asks once more; `--dry-run` only reports). Step 3 uses the GTDB species representatives
(`gtdb_representative == t`): 1,812 genomes, one per species.

The GTDB family Mycobacteriaceae is wider than the genus: 850 of the representatives are *Mycobacterium* (which in GTDB
includes *Mycobacteroides*, *Mycolicibacterium* and the other split genera); the rest are *Corynebacterium*, *Nocardia*,
*Gordonia*, *Rhodococcus* and relatives. `coverage.tsv` counts products for the genus and for the family.

## Reference sequences

A product is reported when both primers match with at most 3 substitutions each (IUPAC codes of the primer honoured,
no gaps) and it is at most 3 kb long. Per genome and pair only the products with the fewest mismatches are kept: weaker
ones are other loci (the hsp65 primers also hit the groEL1 paralog with 3 + 3 mismatches). Identical copies of a genome
are merged (`copies`). The sequence runs from the forward to the reverse primer and holds the genome's own bases under
the primers, not the primer sequence. The mismatches of each record are in the FASTA header and in `amplicons.tsv`;
filter on them when a strict reference set is needed.

## Extra genomes

GTDB merges some clinically distinct taxa into one species and then has no representative for them: in r232
*M. ulcerans* sits in *M. marinum*. `extra_genomes.tsv` (curated: accession, label, reason) adds such genomes from the
shared directory to steps 3 and 4; their label is used as the species and their gene records carry `extra=t`.

## Full-length genes

`genes.tsv` (curated) names per gene its anchor amplicons (tried in order), the feature type (CDS / rRNA), a name pattern
(sanity flag `name_ok` only), the H37Rv locus and, for groEL2, the paralog groEL1. Step 4 locates the gene through the
product of the anchor, never through annotation names: the GFF3 feature of that type overlapping the product is
exported in full, 5'->3' in gene orientation. Without GFF or overlapping feature (~900 genomes have no GFF) the H37Rv
gene is aligned to the contig around the product (`source=homology`). A record always contains the whole amplicon; the
gene is at `gene_start..gene_end` (1-based) inside it. Every copy of a multi-copy gene (rRNA operons) with a distinct
sequence is a record. `paralog_flag=t` marks a record closer to the H37Rv paralog than to the gene.

Header contract (words 2-3 are the GTDB species, as in the amplicon files):

```
>GCF_002102175.1|groEL2|1 Mycobacterium gastri [gene=groEL2 length=1623 gene_start=1 gene_end=1623 gene_bp=1623 source=gff anchor=hsp65_TB11_TB12w copies=1 type_strain=t ncbi_organism="Mycobacterium gastri"]
```

The record id `accession|gene|n` differs from the amplicon ids (`accession|amplicon|n`) of the same genome.

## SNP tables

`snp_pairs.tsv` (curated) lists species pairs: gene, species A / B, the window amplicon (the region that is sequenced;
empty = whole gene), use (`diagnostic` / `research`), species whose allele is reported, and optionally legacy GenBank
records with positions (0-based) used before. Step 5 aligns the representative gene records of A and B, takes the
substitutions inside the window (a pair with an indel in the window is reported but not used), and counts the alleles in
the gene of **every** genome of both species on disk. `diagnostic=t` when both species carry their allele in at least
`--min-fixed` (0.99) of the genomes with a call. Positions are given 0-based in the A and B records and as H37Rv gene
nucleotide / codon. `continuity.tsv` places the legacy positions in the new coordinates (the legacy records are fetched
once from NCBI into `snps/legacy/`).

## Layout

```
cosmic_seq/       handbook.py (text of .docx / .doc), insilico.py (primer search, products), tables.py (tables, step-03 cache),
                  genomes.py (representatives, extra genomes, GTDB metadata), genes.py (GFF3, gene extraction, alignment),
                  kansasii_paths.py (vendored)
scripts/          01_gtdb_mycobacteriaceae_genomes.sh, 02_handbook_primers.py, 03_reference_amplicons.py,
                  04_reference_genes.py, 05_snp_pairs.py
primers.tsv  amplicons.tsv  genes.tsv  extra_genomes.tsv  snp_pairs.tsv   curated tables
docs/PRIMERS.md   generated
submit_references.sbatch     steps 2 to 5
```

Paths are constants of `kansasii_paths` (`HANDBOOK_DIR`, `NCBI_DATASET`, `OUT_COSMIC`) and flags of the scripts.

## Run

```bash
conda env create -f environment.yml          # env_cosmic, once
conda activate env_cosmic
bash scripts/01_gtdb_mycobacteriaceae_genomes.sh [--dry-run]
python scripts/02_handbook_primers.py
sbatch submit_references.sbatch              # 02 -> 05 on 16 cores; log in output/cosmic-seq/logs/
python scripts/03_reference_amplicons.py --accession GCF_000195955.2   # one genome, about 20 s
python scripts/04_reference_genes.py --accession GCF_000195955.2 --out <scratch dir>   # test without touching results
```

Steps 3, 4 and 5 keep finished genomes in `per_genome.jsonl`; a rerun only processes new genomes. Step 3 redoes
everything when a primer, a pairing or a search limit changes; step 4 redoes a genome when `genes.tsv`, the H37Rv
references or its anchor products change. A run with `--accession` rewrites the tables for those genomes only: point
`--out` elsewhere for tests. Type check: `~/data/conda/envs/env_typecheck/bin/mypy`.
