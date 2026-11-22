#!/usr/bin/env bash
# 02 + 03 + 05: in-silico PCR of all primer pairs on the GTDB Mycobacteriaceae representatives -> output/cosmic-seq/results/references/
#SBATCH --job-name=cosmic-references
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
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
python scripts/05_export_cosmic_files.py      # bundle for cosmic-files (needs the species groups of 04)
