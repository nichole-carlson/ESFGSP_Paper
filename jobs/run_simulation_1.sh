#!/bin/bash
#SBATCH --job-name=sim1
#SBATCH --output=logs/sim1_%A_%a.log
#SBATCH --ntasks=1                       # Number of tasks (1 node)
#SBATCH --cpus-per-task=1                # Number of CPU cores per task
#SBATCH --mem=1G                         # Memory per node
#SBATCH --partition=amilan               # Partition for array tasks
#SBATCH --time=00:20:00                  # Max runtime
#SBATCH --array=1-500
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=siyang.ren@cuanschutz.edu

module load R

OUT_DIR=/scratch/alpine/sren@xsede.org/esfgsp/sim1_iters

SIM_ID=${SLURM_ARRAY_TASK_ID}
SEED=$((10000 + SIM_ID))

Rscript --no-save scripts/simulation_1.R \
  --sim_id "$SIM_ID" \
  --effect 0.1 \
  --n_sample 1000 \
  --out_dir "$OUT_DIR" \
  --seed "$SEED"
