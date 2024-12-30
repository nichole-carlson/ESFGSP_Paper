#!/bin/bash
#SBATCH --job-name="model"
#SBATCH --output=logs/model1a_%j.out
#SBATCH --ntasks=1                       # Number of tasks (1 node)
#SBATCH --cpus-per-task=8                # Number of CPU cores per task
#SBATCH --mem=120G                       # Memory per node
#SBATCH --partition=amilan               # Partition for array tasks
#SBATCH --time=24:00:00                  # Max runtime
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=siyang.ren@cuanschutz.edu

# Load modules
module load R

# Run R script to simulate data
Rscript --no-save scripts/02_fit_models.R
