# WRAP MODEL FITTING WITH LEAVE FUTURE OUT CV

library(parallel)
library(caret)
library(parallelly)
library(blockCV)
print(paste("Cores:", detectCores()))
print(paste("Cores:", availableCores()))

suppressMessages(source("code/setup.R"))
suppressMessages(source("code/build_design_matrix.R"))
suppressMessages(source("code/betabinomial_p_rho.R"))
suppressMessages(source("code/wrap_fit.R"))

args <- commandArgs(trailingOnly = TRUE)
marker <- args[1]  # of "k13_marcse", "mdr86", "mdr184", "mdr1246", "crt76" etc
seed <- as.numeric(args[2])
print(paste0("Marker: ", marker))
print(paste0("Seed: ", seed))

set.seed(seed)

out_dir <- paste0(marker, "/bb_gne/lfo/")

in_dat <- data_path_lookup[[marker]]


print(paste0("Reading in from: ", in_dat))
print("Enforcing min year for surveyor data - 2000")
print("Fitting betabinom")

mut_data <- setup_mut_data(in_dat, 
                           min_year = MIN_YEAR, 
                           buffer = BUFFER)
write_rds(mut_data, paste0("output/", out_dir, "mut_data.rds"))

NFOLD <- 10



# separate out folds by time
folds <- lapply(1: (NFOLD + 1), function(x){
  cutoff <- max(mut_data$year) - NFOLD + x - 1
  which(mut_data$year <= cutoff)
})

names(folds) <- paste0("Fold", 1: (NFOLD + 1))

write_rds(folds, paste0("output/", out_dir, "cv_folds_lfo.rds"))

# system.time(mclapply(1:NFOLD, function(x){
#   fit_betabinom(mut_data = mut_data,
#                 covariates = covariates,
#                 pfpr_years = pfpr_years,
#                 out_dir = out_dir,
#                 fold = x,
#                 folds = folds,
#                 lfo = TRUE,
#                 nchains = 6,
#                 warmup = 5000,
#                 nsamples = 30000)
# }, mc.cores = NFOLD))

system.time(mclapply(1:NFOLD, function(x){
  fit_betabinom(mut_data = mut_data,
                covariates = covariates,
                pfpr_years = pfpr_years,
                out_dir = out_dir,
                fold = x,
                folds = folds,
                lfo = TRUE,
                nchains = 6,
                warmup = 5000,
                nsamples = 30000)
}, mc.cores = NFOLD))

