#!/usr/bin/env bash
# 01: GTDB Mycobacteriaceae genomes -> shared $GTDB_DIR (data/gtdb_genomes/Mycobacteriaceae), only what is missing.
#
# Same logic as immensekansasii/scripts/01_gtdb_mycobacteria_genomes.sh (GTDB metadata -> f__Mycobacteriaceae
# accessions -> NCBI `datasets` in chunks), without the kansasii-specific lists and with a stricter
# "already present" test. Both scripts fill the same directory, so whichever runs second downloads nothing.
#
# Lists are written to $OUT_COSMIC/results/gtdb<release>/:
#   mycobacteriaceae_rows_metadata.tsv       all GTDB metadata columns of the family
#   mycobacteriaceae_selected_columns.tsv    accession, representative, taxonomy, type-strain and NCBI columns
#   mycobacteriaceae_accessions.txt          bare GCA_/GCF_ accessions
#   mycobacteriaceae_representatives.tsv     accession, genus, species of the GTDB species representatives (input of 03)
#   unavailable_accessions.txt               requested at NCBI and not delivered; skipped on later runs
#   not_downloaded_genomes.tsv               selected columns of every accession without a genome
#
# Usage: 01_gtdb_mycobacteriaceae_genomes.sh [--retry-unavailable] [--dry-run]
#   --retry-unavailable   request the accessions of unavailable_accessions.txt again
#   --dry-run             write the lists and report what would be downloaded, download nothing
# Needs NCBI `datasets` only when something is missing (conda activate env_cosmic).

set -euo pipefail
export LC_ALL=C  # one collation for sort and comm
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"

RETRY_UNAVAILABLE=0
DRY_RUN=0
for arg in "$@"; do
  case "$arg" in
    --retry-unavailable) RETRY_UNAVAILABLE=1 ;;
    --dry-run) DRY_RUN=1 ;;
    *) echo "unknown argument: $arg" >&2; exit 2 ;;
  esac
done

OUT_DIR=$OUT_COSMIC/results/gtdb$GTDB_RELEASE
CHUNK_SIZE=500
MAX_RETRIES=3
mkdir -p "$GTDB_METADATA_DIR" "$OUT_DIR" "$NCBI_DATASET"

# ---------------------------------------------------------------------------------------------- GTDB metadata
METADATA="$GTDB_METADATA_DIR/bac120_metadata_r$GTDB_RELEASE.tsv"
if [ ! -s "$METADATA" ]; then
  echo "-- downloading GTDB r$GTDB_RELEASE metadata --"
  url="https://data.gtdb.aau.ecogenomic.org/releases/release$GTDB_RELEASE/$GTDB_RELEASE.0/bac120_metadata_r$GTDB_RELEASE.tsv.gz"
  wget -q -O "$METADATA.gz.part" "$url"
  gunzip -c "$METADATA.gz.part" > "$METADATA.part"
  mv "$METADATA.part" "$METADATA"
  rm -f "$METADATA.gz.part"
fi

cd "$OUT_DIR"

# ------------------------------------------------------------------------------------------------------ lists
(head -1 "$METADATA"; grep -F 'f__Mycobacteriaceae' "$METADATA") > mycobacteriaceae_rows_metadata.tsv

column_number() {  # 1-based position of a column of the metadata header
  local n
  n=$(head -1 "$METADATA" | tr '\t' '\n' | grep -n -x "$1" | cut -d: -f1)
  [ -n "$n" ] || { echo "ERROR: column '$1' not in $METADATA" >&2; exit 1; }
  echo "$n"
}
ACCESSION=$(column_number accession)
REPR_COL=$(column_number gtdb_genome_representative)
GTDB_REP_COL=$(column_number gtdb_representative)   # t or f
TAX_COL=$(column_number gtdb_taxonomy)
GTDBTYPE_COL=$(column_number gtdb_type_designation_ncbi_taxa)
ASSEMNAME_COL=$(column_number ncbi_assembly_name)
BIOSAMPLE_COL=$(column_number ncbi_biosample)
NCBITYPE_COL=$(column_number ncbi_type_material_designation)

awk -F'\t' -v OFS='\t' \
  -v acc="$ACCESSION" -v repr="$REPR_COL" -v gtdbrep="$GTDB_REP_COL" \
  -v tax="$TAX_COL" -v gtype="$GTDBTYPE_COL" -v an="$ASSEMNAME_COL" \
  -v bs="$BIOSAMPLE_COL" -v ntype="$NCBITYPE_COL" \
  '{print $acc, $repr, $gtdbrep, $tax, $gtype, $an, $bs, $ntype}' \
  mycobacteriaceae_rows_metadata.tsv > mycobacteriaceae_selected_columns.tsv

# GTDB accessions carry a 3-character source flag ("RS_"/"GB_") that NCBI does not accept: strip it
tail -n +2 mycobacteriaceae_rows_metadata.tsv | cut -f"$ACCESSION" | cut -c4- | sort -u > mycobacteriaceae_accessions.txt
[ -s mycobacteriaceae_accessions.txt ] || { echo "ERROR: no f__Mycobacteriaceae rows in $METADATA" >&2; exit 1; }

# GTDB species representatives of the whole family (gtdb_representative == t), one genome per species.
# The representative is chosen by assembly quality and is often not the type strain.
awk -F'\t' -v OFS='\t' -v acc="$ACCESSION" -v tax="$TAX_COL" -v gtdbrep="$GTDB_REP_COL" '
  NR == 1 {print "accession", "genus", "species"; next}
  $gtdbrep == "t" {
    n = split($tax, rank, ";"); genus = ""; species = ""
    for (i = 1; i <= n; i++) {
      if (rank[i] ~ /^g__/) genus = substr(rank[i], 4)
      if (rank[i] ~ /^s__/) species = substr(rank[i], 4)
    }
    print substr($acc, 4), genus, species
  }' mycobacteriaceae_rows_metadata.tsv > mycobacteriaceae_representatives.tsv

# --------------------------------------------------------------------------------------------------- download
present_accessions() {  # accessions whose directory holds a non-empty genome FASTA (one directory scan)
  find "$NCBI_DATASET" -mindepth 2 -maxdepth 2 -name '*_genomic.fna' -size +0 -printf '%h\n' \
    | xargs -r -n 1000 basename -a | sort -u
}

fetch_chunk() {  # fetch_chunk ACCESSION_FILE: one `datasets` call, unpacked into $GTDB_DIR
  # Only the genome directories are unpacked: README.md, md5sum.txt and the package-level reports of the
  # archive describe one request and would overwrite those of the previous one.
  local accfile=$1 zip="$1.zip" status=0
  if datasets download genome accession --inputfile "$accfile" --include gff3,genome --no-progressbar --filename "$zip"; then
    unzip -o -q "$zip" 'ncbi_dataset/data/GC*' -d "$GTDB_DIR" || status=$?
    [ "$status" -eq 11 ] && status=0   # 11 = archive holds no genome directory (records without sequence)
  else
    status=1
  fi
  rm -f "$zip"
  return "$status"
}

download_missing() {
  local chunkdir="$GTDB_DIR/.accession_chunks_cosmic"
  local present="$chunkdir/present.txt" missing="$chunkdir/missing.txt" todo="$chunkdir/todo.txt"
  mkdir -p "$chunkdir"
  touch unavailable_accessions.txt

  present_accessions > "$present"
  comm -23 mycobacteriaceae_accessions.txt "$present" > "$missing"
  if [ "$RETRY_UNAVAILABLE" -eq 1 ]; then
    cp "$missing" "$todo"
  else
    comm -23 "$missing" <(sort -u unavailable_accessions.txt) > "$todo"
  fi

  local n_all n_missing n_todo
  n_all=$(wc -l < mycobacteriaceae_accessions.txt); n_missing=$(wc -l < "$missing"); n_todo=$(wc -l < "$todo")
  echo "-- Mycobacteriaceae: $n_all accession(s), $((n_all - n_missing)) present, $n_missing missing ($((n_missing - n_todo)) known unavailable at NCBI) --"
  if [ "$n_todo" -eq 0 ]; then
    echo "-- nothing to download --"
    return 0
  fi
  if [ "$DRY_RUN" -eq 1 ]; then
    echo "-- dry run: would request $n_todo accession(s), listed in $todo --"
    return 0
  fi
  command -v datasets > /dev/null || { echo "ERROR: NCBI 'datasets' CLI not found (conda activate env_cosmic)" >&2; return 1; }
  command -v unzip > /dev/null || { echo "ERROR: unzip not found" >&2; return 1; }

  rm -f "$chunkdir"/chunk_*
  echo "-- requesting $n_todo accession(s) in chunks of $CHUNK_SIZE --"
  split -d -a 4 -l "$CHUNK_SIZE" "$todo" "$chunkdir/chunk_"

  local chunk attempt ok acc
  for chunk in "$chunkdir"/chunk_*; do
    ok=0
    for attempt in $(seq 1 "$MAX_RETRIES"); do
      echo "-- chunk $(basename "$chunk"): attempt $attempt/$MAX_RETRIES ($(wc -l < "$chunk") accession(s)) --"
      if fetch_chunk "$chunk"; then ok=1; break; fi
      echo "chunk $(basename "$chunk") failed on attempt $attempt" >&2
    done
    if [ "$ok" -ne 1 ]; then
      # `datasets` fails for a whole request when none of its accessions exists any more (suppressed or
      # replaced records). Ask for each accession on its own so that one dead record cannot block the rest.
      echo "-- chunk $(basename "$chunk"): requesting its accessions one by one --" >&2
      while read -r acc; do
        echo "$acc" > "$chunkdir/single"
        fetch_chunk "$chunkdir/single" || echo "   $acc: not delivered" >&2
      done < "$chunk"
      rm -f "$chunkdir/single"
    fi
  done

  # Requested and still without a genome: remembered, so that the next run does not ask again.
  # NCBI must be reachable for this to mean "unavailable": if it is not, nothing is recorded.
  present_accessions > "$present"
  local still="$chunkdir/still_missing.txt"
  comm -23 "$todo" "$present" > "$still"
  if [ -s "$still" ]; then
    if datasets summary genome taxon 'Mycobacterium tuberculosis' --limit 1 > /dev/null 2>&1; then
      sort -u "$still" unavailable_accessions.txt > "$chunkdir/unavailable.new"
      mv "$chunkdir/unavailable.new" unavailable_accessions.txt
      echo "-- $(wc -l < "$still") requested accession(s) not delivered by NCBI, added to unavailable_accessions.txt --"
    else
      echo "ERROR: NCBI not reachable; $(wc -l < "$still") accession(s) still missing, run again" >&2
      return 1
    fi
  fi
  # an accession that arrived after all is no longer unavailable
  comm -23 <(sort -u unavailable_accessions.txt) "$present" > "$chunkdir/unavailable.new"
  mv "$chunkdir/unavailable.new" unavailable_accessions.txt
  rm -f "$chunkdir"/chunk_*
}

download_missing

if [ "$DRY_RUN" -eq 0 ]; then
  # every accession without a genome (suppressed or withdrawn at NCBI), same columns as the selected table
  present_accessions > "$GTDB_DIR/.accession_chunks_cosmic/present.txt"
  awk -F'\t' -v OFS='\t' 'FILENAME == ARGV[1] {have[$1]; next} FNR == 1 || !(substr($1, 4) in have)' \
    "$GTDB_DIR/.accession_chunks_cosmic/present.txt" mycobacteriaceae_selected_columns.tsv > not_downloaded_genomes.tsv
  echo "-- $(($(wc -l < not_downloaded_genomes.tsv) - 1)) genome(s) not available, see $OUT_DIR/not_downloaded_genomes.tsv --"
  n_rep=$(($(wc -l < mycobacteriaceae_representatives.tsv) - 1))
  n_rep_missing=$(awk -F'\t' 'FILENAME == ARGV[1] {have[$1]; next} FNR > 1 && !($1 in have)' \
    "$GTDB_DIR/.accession_chunks_cosmic/present.txt" mycobacteriaceae_representatives.tsv | wc -l)
  echo "-- $n_rep species representative(s), $n_rep_missing of them without a genome --"

  # provenance (see kansasii-docs/docs/PROVENANCE.md)
  bash "$MANIFEST_SCRIPT" "$OUT_DIR" "$REPOS/cosmic-seq" scripts/01_gtdb_mycobacteriaceae_genomes.sh "$METADATA"
fi
