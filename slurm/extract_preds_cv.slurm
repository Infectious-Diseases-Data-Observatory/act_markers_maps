#!/bin/bash
#SBATCH --job-name=extr_pre
#SBATCH --output=logs/%x_%A_.out
#SBATCH --error=logs/%x_%A_.err
#SBATCH --time=1:00:00           # Max runtime (HH:MM:SS)
#SBATCH --mem=5G                  # Memory per node
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --cpus-per-task=1         # Number of CPU cores per task

#SBATCH --mail-type=BEGIN,END,FAIL  # Send email when job begins, ends, or fails
#SBATCH --mail-user=lucy.harrison@wwarn.org  # Replace with your email address

# should take under 1 h

module load Miniforge3/24.1.2-0
eval "$(conda shell.bash hook)"
conda activate greta-env

cd ../

pwd

# provide marker and model ...
Rscript code/extract_preds_cv.R k13_marcse bb_gne
Rscript code/extract_preds_cv.R mdr86 bb_gne
Rscript code/extract_preds_cv.R mdr184 bb_gne
Rscript code/extract_preds_cv.R mdr1246 bb_gne
Rscript code/extract_preds_cv.R crt76 bb_gne
Rscript code/extract_preds_cv.R k13snp_R561H bb_gne
Rscript code/extract_preds_cv.R k13snp_R622I bb_gne
Rscript code/extract_preds_cv.R k13snp_A675V bb_gne
Rscript code/extract_preds_cv.R k13snp_P441L bb_gne
Rscript code/extract_preds_cv.R k13snp_C469Y bb_gne




echo "Job done!"
