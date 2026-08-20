#!/bin/bash
#SBATCH --job-name=inf_cv_lfo
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=4:00:00           # Max runtime (HH:MM:SS)
#SBATCH --mem=50000                  # Memory per node
#SBATCH --array=0-2
#SBATCH --cpus-per-task=30         # Number of CPU cores per task

#SBATCH --mail-type=BEGIN,END,FAIL  # Send email when job begins, ends, or fails
#SBATCH --mail-user=lucy.harrison@wwarn.org  # Replace with your email address

module load Miniforge3/24.1.2-0
eval "$(conda shell.bash hook)"
conda activate greta-env

cd ../

pwd

# mop up some dodgy LFO models for partner drug datasets ...

DATASETS=(
  mdr86
  mdr184
  mdr184
)

SEEDS=(
  124
  124
  124
)

# a number between 1 and 10
FOLDS=(
  6
  4
  5
)

# Select dataset for this array task
DATASET=${DATASETS[$SLURM_ARRAY_TASK_ID]}

SEED=${SEEDS[$SLURM_ARRAY_TASK_ID]}

FOLD=${FOLDS[$SLURM_ARRAY_TASK_ID]}

# Run R script
Rscript code/fit_betabinom_cv_lfo.R "$DATASET" "$SEED" "$FOLD"

echo "Job done!"
