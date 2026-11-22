# cosmic-seq

Reference sequences for the mycobacterial Sanger assays of the lab: the primers documented in the lab handbooks
(`data/Arbeitshandbücher`) are run in silico on the GTDB species representatives of the family Mycobacteriaceae, and
each product is kept as a reference for the reads of that assay.

Five steps:

| Step | Script | Does | Output |
|---|---|---|---|
| 1 | `scripts/01_gtdb_mycobacteriaceae_genomes.sh` | GTDB r232 metadata -> all `f__Mycobacteriaceae` accessions -> NCBI `datasets`, **only genomes that are not on disk** | genomes in the shared `data/gtdb_genomes/Mycobacteriaceae/ncbi_dataset/data/`; lists in `output/cosmic-seq/results/gtdb232/` |
| 2 | `scripts/02_handbook_primers.py` | checks `primers.tsv` against the protocol files, lists other oligos of the mycobacterial protocols, runs every pair on H37Rv | `docs/PRIMERS.md`; `output/cosmic-seq/results/primers/` |
| 3 | `scripts/03_reference_amplicons.py` | in-silico PCR of every primer pair on every species representative | `output/cosmic-seq/results/references/` (`fasta/<amplicon>.fasta`, `amplicons.tsv`, `coverage.tsv`) |
| 4 | `scripts/04_ncbi_species_groups.py` | species groups (complexes): NCBI taxonomy plus the curated `species_groups.tsv`, one group per species | `output/cosmic-seq/results/species_groups/` (`species_groups.tsv`, `groups.tsv`, `raw/`) |
| 5 | `scripts/05_export_cosmic_files.py` | bundle for cosmic-files: the amplicon FASTAs of 3 and `species_groups.tsv` of 4 | `output/cosmic-seq/results/cosmic_files/` |

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

## Species groups

Step 4 asks the NCBI Datasets taxonomy API for every `SPECIES_GROUP` / `SPECIES_SUBGROUP` below Mycobacteriaceae (taxid
1762; found, not listed by hand) and for the species of each (`filtered_subtree?rank_limits=SPECIES`). In November 2026
these are the M. tuberculosis complex, the M. avium complex, the M. simiae complex and the M. ulcerans group; NCBI has no
kansasii, abscessus, fortuitum or terrae complex. NCBI "species" include placeholders (`Mycobacterium sp. ...`, `cf.`,
quoted isolate names, mostly without assemblies); they are counted in `groups.tsv` but left out of `species_groups.tsv`.

The curated `species_groups.tsv` of the repo (`group, label, species, note`) adds the groups NCBI lacks (M. kansasii
complex, M. chelonae / M. abscessus complex, with both *Mycobacterium* and *Mycobacteroides* spellings) and gives each
group its display label; a row without `species` only sets the label of an NCBI group. A curated species that is also in
an NCBI group stops the script unless its note says `override`.

The output `species_groups.tsv` (`group, label, species, source, group_taxid, species_taxid, ncbi_assemblies`) holds
one row per formal species and one group per species: a lookup "species name -> group" for cosmic-files, which names a
16S call by its group because 16S does not separate the members of a complex. GTDB is not used: it has no rank for
complexes. The script needs network (about 5 requests, about 20 s) and no cluster data, so it also runs on Windows;
`--from-raw` parses the saved answers in `raw/` again. Without bash or the manifest script it skips `MANIFEST.txt`.

## Export for cosmic-files

cosmic-files embeds the references at compile time from its `res/sequences/references/`. Step 5 writes exactly that
directory, flat: `<amplicon>.fasta`, `species_groups.tsv`, `MANIFEST.txt`; FASTA files that step 3 no longer writes are
removed. It checks the FASTA headers and that no species is in two groups. On the desktop, cosmic-files copies the
bundle with `python scripts/sync_cosmic_seq.py --remote <host>:<results>/cosmic_files` (scp) or `--from <dir>`.
Step 5 needs no genome data and also runs on Windows (`--results <copy of results>`).

## Layout

```
cosmic_seq/       handbook.py (text of .docx / .doc), insilico.py (primer search, products), tables.py, kansasii_paths.py (vendored)
scripts/          01_gtdb_mycobacteriaceae_genomes.sh, 02_handbook_primers.py, 03_reference_amplicons.py, 04_ncbi_species_groups.py,
                  05_export_cosmic_files.py
primers.tsv  amplicons.tsv  species_groups.tsv   curated tables
docs/PRIMERS.md   generated
submit_references.sbatch     steps 2, 3 and 5
```

Paths are constants of `kansasii_paths` (`HANDBOOK_DIR`, `NCBI_DATASET`, `OUT_COSMIC`) and flags of the scripts.

## Run

```bash
conda env create -f environment.yml          # env_cosmic, once
conda activate env_cosmic
bash scripts/01_gtdb_mycobacteriaceae_genomes.sh [--dry-run]
python scripts/02_handbook_primers.py
sbatch submit_references.sbatch              # 02 + 03 + 05 on 16 cores; log in output/cosmic-seq/logs/
python scripts/03_reference_amplicons.py --accession GCF_000195955.2   # one genome, about 20 s
python scripts/04_ncbi_species_groups.py     # species groups, needs network
python scripts/05_export_cosmic_files.py     # bundle for cosmic-files (also run by the sbatch)
```

Step 3 keeps finished genomes in `per_genome.jsonl`; a rerun only processes new genomes, and everything is redone when a
primer, a pairing or a search limit changes. Type check: `~/data/conda/envs/env_typecheck/bin/mypy`.
