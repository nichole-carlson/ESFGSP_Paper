#!/bin/bash
#SBATCH --job-name=s2_com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --partition=amilan
#SBATCH --time=00:30:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=siyang.ren@cuanschutz.edu
#SBATCH --output=logs/sim2_combine_%j.log

module load R

DATA_DRIVE="/scratch/alpine/sren@xsede.org/esfgsp"

RESULT_DIR="${DATA_DRIVE}/sim2_iters"
OUT_DIR="${DATA_DRIVE}"

Rscript --no-save scripts/combine_sim_results.R \
        --indir $RESULT_DIR \
        --outdir $OUT_DIR \
        --sim_id "sim2"
