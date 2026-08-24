# extract predictions for held-out runs
start <- Sys.time()

suppressMessages(source("code/setup.R"))
suppressMessages(source("code/validation_funcs.R"))

args <- commandArgs(trailingOnly = TRUE)
marker <- args[1]
mod <- args[2]
message(paste0("Marker: ", marker))
message(paste0("Model: ", mod))

if (!str_detect(mod, "lfo")){
    folds <- read_rds(paste0("output/", marker, "/", mod, "/cv_folds_spat.rds"))
    pred_path = paste0("output/", marker, "/", mod, "/cv_preds/")
    lfo  <- FALSE
} else {
    # mod includes "/lfo"
    folds <- read_rds(paste0("output/", marker, "/", mod, "/cv_folds_lfo.rds"))
    pred_path = paste0("output/", marker, "/", mod, "/")
    lfo <- TRUE
}


dat <- extract_preds_cv(data_path = data_path_lookup[[marker]],
                 pred_path = pred_path,
                 folds = folds, 
                 lfo = lfo,
                 in_buffer = BUFFER,
                 out_buffer = 50000) 
# using large buffer as increased AGG_FACTOR for held-out models
# this completely removed some islands from the raster and chopped up some of the internal 
# waterbodies/holes in the transmission mask
# ultimately we have explicit spatial covariance so it ain't a huge problem if we need to reland 
# some points

write.csv(dat, paste0("output/", marker, "/", mod, "/mut_dat_cv_preds_extracted.csv"),
          row.names = FALSE)

end <- Sys.time()

print(end - start)