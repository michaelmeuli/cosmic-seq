#!/usr/bin/env bash
# 02 -> 05: primer check, in-silico PCR (references/), full-length genes (genes/), diagnostic SNP tables (snps/)
# on the GTDB Mycobacteriaceae representatives + extra_genomes.tsv -> output/cosmic-seq/results/
#SBATCH --job-name=cosmic-references
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/cosmic-seq/logs/%x-%j.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/cosmic-seq
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${COSMIC_ENV:-env_cosmic}"

cd "$REPO"
# single-threaded numpy per worker: the workers are the parallelism
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python scripts/02_handbook_primers.py          # stops here if a primer is not in its protocol
python scripts/03_reference_amplicons.py --jobs "${SLURM_CPUS_PER_TASK:-1}" "$@"
python scripts/04_reference_genes.py --jobs "${SLURM_CPUS_PER_TASK:-1}"
python scripts/05_snp_pairs.py --jobs "${SLURM_CPUS_PER_TASK:-1}"
