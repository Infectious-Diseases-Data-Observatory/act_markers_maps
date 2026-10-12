set.seed(0784)

# TODO
# filter on end.year ..
# filter on 

suppressPackageStartupMessages({
  library(greta)
  library(greta.gp)
  library(terra)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(magrittr)
  library(tidyverse)
  library(sf) # not sure why I'm adding this in October 2025
  library(rnaturalearth)
  library(rnaturalearthdata)})
world <- ne_countries(scale="medium", returnclass = "sf")

MIN_YEAR <- 2000
BUFFER <- 5000

theme_set(theme_bw())

# here is some run up:
nice_name_lookup_main <- list("k13_marcse" = "Kelch 13",
                              "crt76" = "Pfcrt-K76T",
                              "mdr86" = "Pfmdr1-N86Y",
                              "mdr184" = "Pfmdr1-Y184F",
                              "mdr1246" = "Pfmdr1-D1246Y")

nice_name_lookup_all <- list("k13_marcse" = "Kelch 13 aggregate",
                             "crt76" = "Pfcrt-K76T",
                             "mdr86" = "Pfmdr1-N86Y",
                             "mdr184" = "Pfmdr1-Y184F",
                             "mdr1246" = "Pfmdr1-D1246Y",
                             "k13snp_A675V" = "Kelch 13 A675V",
                             "k13snp_C469Y" = "Kelch 13 C469Y",
                             "k13snp_P441L" = "Kelch 13 P441L",
                             "k13snp_R561H" = "Kelch 13 R561H",
                             "k13snp_R622I" = "Kelch 13 R622I")

data_path_lookup <- list("k13_marcse" = "data/clean/moldm_marcse_k13_nomarker.csv",
                         "crt76" = "data/clean/moldm_crt76.csv",
                         "mdr86" = "data/clean/pfmdr_single_mdr86.csv",
                         "mdr184" = "data/clean/pfmdr_single_mdr184.csv",
                         "mdr1246" = "data/clean/pfmdr_single_mdr1246.csv",
                         "k13snp_A675V" = "data/clean/moldm_marcse_k13snp_A675V.csv",
                         "k13snp_C469Y" = "data/clean/moldm_marcse_k13snp_C469Y.csv",
                         "k13snp_P441L" = "data/clean/moldm_marcse_k13snp_P441L.csv",
                         "k13snp_R561H" = "data/clean/moldm_marcse_k13snp_R561H.csv",
                         "k13snp_R622I" = "data/clean/moldm_marcse_k13snp_R622I.csv")

bb_paths <- lapply(names(nice_name_lookup_all),
                   function(marker){paste0("output/", marker, "/bb_gne/")})
names(bb_paths) <- names(nice_name_lookup_all)


# country shps for plotting/masking
afr <- world %>%
  filter(continent == "Africa") %>%
  vect() %>%
  crop(ext(-21, 63, -35, 37)) %>%
  st_as_sf()

# make some edits to rnaturalearth subregion labels
# Zambia to Southern Africa, Angola to Southern Africa
# Sudan to Eastern Africa, "Middle Africa" to "Central Africa"
afr <- afr %>%
  mutate(subregion = case_when(name_en == "Zambia" ~ "Southern Africa",
                               name_en == "Angola" ~ "Southern Africa",
                               name_en == "Sudan" ~ "North-eastern Africa",
                               name_en == "Ethiopia" ~ "North-eastern Africa",
                               name_en == "South Sudan" ~ "North-eastern Africa",
                               name_en == "Eritrea" ~ "North-eastern Africa",
                               subregion == "Middle Africa" ~ "Central Africa",
                               TRUE ~ subregion))
# (this is for colouring points in ribbon plots)

afr_extent <- (matrix(c(-21, -34.9, 63, 37.4), nrow = 2, ncol = 2, 
                      dimnames = list(c("x", "y"), c("min", "max"))))

# For example: 
pfpr_years <- 2000:2024
years <- as.character(pfpr_years)

# see "retrieve_pfpr.R"
pfpr <- rast("data/pfpr_rasters_afr_2025.tif") %>%
  scale()
names(pfpr) <- paste0("pfpr_", years)

setup_mut_data <- function(path, min_year = NULL, buffer = 0){
  
  # read in and format data for a single response (k13, crt76, mdr1-86, ...)
  read.csv(path) %>% 
    rename(x = Longitude,
           y = Latitude,
           present = Present,
           tested = Tested) %>%
    filter(if (is.null(min_year)) TRUE else year >= min_year) %>%
    filter(x < afr_extent["x","max"] & afr_extent["x","min"] < x &
             y < afr_extent["y","max"] & afr_extent["y","min"] < y) %>%
    dplyr::select(x, y, year, tested, present) %>%
    group_by(x, y, year) %>%
    summarise(tested = sum(tested),
              present = sum(present)) %>%
    ungroup() %>%
    mutate(land = terra::extract(pfpr$pfpr_2000, 
                                 data.frame(x = x, y = y), 
                                 ID = FALSE,
                                 search_radius = buffer)) %>%
    drop_na() %>%
    dplyr::select(-c(land))
}

# (this was for abandoned pfmdr haplotype model:)
setup_multiple_snps <- function(path, min_year = NULL, buffer = 0){
  # read in and format data for multiple responses
  # (extra columns: snp, snp_id; snp included in grouping)
  read.csv(path) %>% 
    rename(x = Longitude,
           y = Latitude,
           present = Present,
           tested = Tested,
           snp = realised) %>%
    # remove floaters outside of rectangular extent
    filter(x < afr_extent["x","max"] & afr_extent["x","min"] < x &
             y < afr_extent["y","max"] & afr_extent["y","min"] < y) %>%
    dplyr::select(x, y, year, tested, present, snp) %>%
    group_by(x, y, year, snp) %>%
    summarise(tested = sum(tested),
              present = sum(present)) %>%
    ungroup() %>%
    mutate(land = terra::extract(pfpr$pfpr_2000, cbind(x, y), 
                                 search_radius = buffer),
           snp_id = as.numeric(as.factor(snp))) %>%
    drop_na() %>%
    dplyr::select(-c(land))
}

# e.g.:
# mut_data <- setup_multiple_snps("../moldm/clean/pfmdr_single_locus.csv")

# for vis:
library(viridisLite)
library(RColorBrewer)
# library(iddoPal) # excising so that this runs on slurm
oranges <- brewer.pal(9, "Oranges")
# blrd <- iddoPal::iddo_palettes$BlGyRd
blrd <- c("#1481E7", "#14B1E7", "#6ECAEF", "#ACE1F6", "grey80",  "#F19BA0", 
          "#E7515A", "#E31B23", "#A7010C")
# I would like to mess around with this and make the darkest blue a bit darker and more grey in the middle
# blrd_plus <- iddoPal::iddo_palettes$BlGyRdPlus

# iddoblue <- iddo_palettes_discrete$iddo[1]
iddoblue <- "#14B1E7"
# taking some bright colours that don't coincide with viridis:
case_pal <- c("#E37210", iddoblue, "#c7047c")









