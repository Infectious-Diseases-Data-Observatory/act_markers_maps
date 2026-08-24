start <- Sys.time()

suppressMessages(source("code/setup.R"))
suppressMessages(source("code/predict_to_raster.R"))

args <- commandArgs(trailingOnly = TRUE)
marker <- args[1]
mod <- args[2]
fold <- args[3]
message(paste0("Marker: ", marker))
message(paste0("Model: ", mod))
message(paste0("Fold: ", fold))

# set this to the location where all the inference outputs are at:
# out_dir <- "output/mdr86/gneiting_sparse/"
if (!str_detect(mod, "lfo")){   
    out_dir <- paste0("output/", marker, "/", mod, "/cv_preds/")
} else {
    out_dir <- paste0("output/", marker, "/", mod)
}

concat_coverages(out_dir, fold)
concat_preds(out_dir, 
             fold = fold,
             medians = TRUE, 
             sds = TRUE, 
             sdscaled = TRUE, 
             ciwidth = TRUE,
             upper = TRUE,
             lower = TRUE)

end <- Sys.time()

print(end - start)