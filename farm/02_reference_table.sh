#!/bin/bash
# ============================================================================
# Step 2 on the UC Davis Farm cluster: the 45,000-simulation reference table,
# as a SLURM array of 150 tasks x 300 simulations (R/02_reference_table.R).
#
# Submit from the repository root:
#   sbatch farm/02_reference_table.sh
# Output: results/reftable_g10_chunk_{1..150}.rds  ->  then run step 3 locally.
#
# Runs on the preemptible low partition with --requeue; each task checkpoints
# in sub-batches, so a preemption only costs the batch in flight. Benchmark:
# ~110 s per simulation single-core, ~70 min per 300-simulation task on 8 cores.
# Requires reCYCLing v1.9.2 in the user R library (farm/install_recycling.sh).
# ============================================================================
#SBATCH --job-name=reftable_g10
#SBATCH --partition=low
#SBATCH --account=publicgrp
#SBATCH --requeue
#SBATCH --array=1-150
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=farm/logs/reftable_g10_%A_%a.out
#SBATCH --error=farm/logs/reftable_g10_%A_%a.out

source /etc/profile.d/modules.sh
module load R/4.4.2
set -euo pipefail
mkdir -p farm/logs results

N_PER_TASK=${N_PER_TASK:-300}
NBATCH=${NBATCH:-10}
CORES=${SLURM_CPUS_PER_TASK:-8}
SEED=${SLURM_ARRAY_TASK_ID:-1}
OUT="results/reftable_g10_chunk_${SEED}.rds"

echo "[task $SEED] N=$N_PER_TASK nbatch=$NBATCH cores=$CORES -> $OUT ($(date); restarts=${SLURM_RESTART_COUNT:-0})"
Rscript R/02_reference_table.R "$N_PER_TASK" "$CORES" "$SEED" "$OUT" "$NBATCH"
echo "[task $SEED] done ($(date))"
