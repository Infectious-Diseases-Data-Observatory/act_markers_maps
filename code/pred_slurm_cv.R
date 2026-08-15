source("code/setup.R")
source("code/build_design_matrix.R")
source("code/predict_to_raster.R")

# for prediction raster:
AGG_FACTOR = 5

# this feels a bit unflashy but I can't keep having separate scripts
args <- commandArgs(trailingOnly = TRUE)
marker <- args[1]
mod <- args[2]
year <- as.numeric(args[3])
fold <- as.numeric(args[4])
message(paste0("Marker: ", marker))
message(paste0("Model: ", mod))
message(paste0("Year: ", year))
message(paste0("Fold: ", fold))

# set this to the location where all the inference outputs are at:
# out_dir <- "output/mdr86/gneiting_sparse/"
out_dir <- paste0("output/", marker, "/", mod, "/")

# message(pfpr_years)
# scaled_years <- scale_years(range(pfpr_years))
# message(scale_years)

# bring in all of the other outputs here too
mut_data <- read_rds(paste0(out_dir, "mut_data.rds"))
folds <- read_rds(paste0(out_dir, "cv_folds_spat.rds"))
train_dat <- mut_data[unlist(folds[-c(fold)]),]
test_dat <- mut_data[unlist(folds[fold]),]

message(nrow(train_dat))
message(names(train_dat))

# this is now how I'm calculating scaled_years during inference - 
# I'm pretty sure we have points in 23 and 24 for all markers ....
message(range(train_dat$year))
scaled_years <- scale_years(range(train_dat$year))
stable_transmission_mask <- rast("data/stable_transmission_mask.grd") %>%
  aggregate(AGG_FACTOR)
random_field <- read_rds(paste0(out_dir, "random_field_", fold, ".rds"))
parameters <- read_rds(paste0(out_dir, "parameters_", fold, ".rds"))
draws <- read_rds(paste0(out_dir, "draws_", fold, ".rds"))

# boo:
if (!year %in% as.numeric(names(scaled_years))){
  scaled_years <- extended_scaled_years(scaled_years, year)
}

set.seed(0748)
preds <- predict_to_ras(covariates,
                        year,
                        draws,
                        parameters,
                        random_field,
                        agg_factor = AGG_FACTOR,
                        scaled_year = scaled_years[[as.character(year)]],
                        coord_cols = c("x_rd", "y_rd", "year_scaled"),
                        design_cols = c("intercept", "year_scaled", "pfpr"),
                        stable_transmission_mask = stable_transmission_mask,
                        coverage = TRUE,
                        test_pts_for_coverage = test_dat)

# perhaps give me a quick plot here?

writeRaster(preds$out, paste0(out_dir, "cv_preds/", year, "_preds_", fold, ".grd"), overwrite = TRUE)
if(!is.null(preds$coverages)){
  write.csv(preds$coverages, paste0(out_dir, "cv_preds/", year, "_coverages_", fold, ".csv"), row.names = FALSE)
}

message(paste0("written to ", out_dir))





