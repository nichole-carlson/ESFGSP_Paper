#!/bin/bash
#SBATCH --job-name=s1_com
#SBATCH --output=logs/sim1_combine_%j.log
#SBATCH --ntasks=1                       # Number of tasks (1 node)
#SBATCH --cpus-per-task=1                # Number of CPU cores per task
#SBATCH --mem=6G                         # Memory per node
#SBATCH --partition=amilan               # Partition for array tasks
#SBATCH --time=00:30:00                  # Max runtime
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=siyang.ren@cuanschutz.edu

module load R

DATA_DRIVE="/scratch/alpine/sren@xsede.org/esfgsp"

RESULT_DIR="${DATA_DRIVE}/sim1_iters"
OUT_DIR="${DATA_DRIVE}"

Rscript --no-save scripts/combine_sim1_results.R \
        --result_dir $RESULT_DIR \ 
        --out_dir $OUT_DIR
