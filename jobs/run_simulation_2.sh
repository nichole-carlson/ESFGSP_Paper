#!/bin/bash
#SBATCH --job-name=sim2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --partition=amilan
#SBATCH --array=1-500
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=siyang.ren@cuanschutz.edu
#SBATCH --output=logs/sim2_%A.log

module load R

OUT_DIR=/scratch/alpine/sren@xsede.org/esfgsp/sim2_iters

SIM_ID=${SLURM_ARRAY_TASK_ID}
SEED=$((10000 + SIM_ID))

Rscript --no-save scripts/simulation_2.R \
        --sim_id $SIM_ID \
        --out_dir "$OUT_DIR" \
        --effect 0.03 \
        --sparsity 0.1 \
        --n_sample 1000 \
        --seed $SEED
