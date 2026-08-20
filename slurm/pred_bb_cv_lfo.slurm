#!/bin/bash
#SBATCH --job-name=pre_bb_cv_lfo
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=3:00:00           # Max runtime (HH:MM:SS)
#SBATCH --mem=60G                  # Memory per node
#SBATCH --array=0-29                 # should be up to nmarkers*nfolds
#SBATCH --cpus-per-task=20         # Number of CPU cores per task

#SBATCH --mail-type=BEGIN,END,FAIL  # Send email when job begins, ends, or fails
#SBATCH --mail-user=lucy.harrison@wwarn.org  # Replace with your email address

module load Miniforge3/24.1.2-0
eval "$(conda shell.bash hook)"
conda activate greta-env

pwd

cd ../

markers=(
  mdr86
  mdr184
  mdr1246
  crt76
  k13_marcse
  # k13snp_R561H
  # k13snp_R622I
  # k13snp_C469Y
  # k13snp_P441L
  # k13snp_A675V
)

nfolds=6

idx=$SLURM_ARRAY_TASK_ID

#fold=$(( idx % ${#markers[@]}  + 1))
fold=$(( idx % nfolds  + 1))

# Select marker for this array task
#marker=${markers[idx / ${#markers[@]}]}
marker=${markers[idx / nfolds]}

echo $idx $marker $fold

make_preds() {
    local MARKER=$1
    local MODEL=$2
    local FOLD=$3
    
    echo "Running task for marker: $MARKER, model: $MODEL, fold: $FOLD"

    # Run the R script with the parameters
    Rscript code/pred_slurm_cv_lfo.R $MARKER $MODEL $FOLD

}

make_preds "$marker" "bb_gne/lfo" "$fold"


echo "Job done! ${SLURM_ARRAY_TASK_ID}"