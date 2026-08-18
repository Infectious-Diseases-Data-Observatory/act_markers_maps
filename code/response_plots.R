
extract_pfpr <- function(df,
                          covs,
                          max_year = 2024,
                          buffer = 0){
  # it's time to throw in the towel chatty g can do this way better than me
  yrs_covs <- str_extract(names(covs), "\\d{4}")
  
  # get predictions for each row in `mut_data`
  df$pfpr <- NA
  yrs_to_extract <- unique(df$year)
  for (yr in yrs_to_extract){
    message(yr)
    idx <- which(df$year == yr)
    
    yr <- min(yr, max_year)
    val <- terra::extract(covs[[paste0("pfpr_", yr)]], 
                          df[idx, c("x", "y")],
                          ID = FALSE, 
                          search_radius = buffer) 
    
    df[idx, "pfpr"] <- val[, paste0("pfpr_", yr)]
  }
  
  df
}

pair_of_rasts <- function(ras, year, pfpr_upper_limit = 2024){
  # give me the right pair of pfpr and prediction rasters for a given year
  c(ras[paste0(year, "_50")],
    ras[paste0("pfpr_", min(as.numeric(year), pfpr_upper_limit))]) %>%
    setNames(c("pred", "pfpr"))
}

#' Run-up to response plots
#' Partitioning as it's a big old object and I don't want to run every time
#' I make up the plot
#'
#' @param ras rast of pfprs and preds
#'
#' @returns df of pfprs and preds
#' @export
#'
#' @examples
response_plots_runup <- function(ras){
  years <- ras %>%
    names() %>%
    str_extract(pattern = "\\d{4}")
  
  lapply(years, function(year){
    pair_of_rasts(ras, year) %>%
      as.data.frame() %>%
      mutate(year = as.numeric(year)) %>%
      drop_na()
  }) %>%
    do.call(what = rbind)
}

#' "Response" plots
#' How does prediction vary over covariate ?
#'
#' @param preds rast
#' @param covar rast
#'
#' @returns
#' @export
#'
#' @examples
response_plot <- function(df, dat = NULL, covar = "pfpr", xax_breaks = 100,
                          pred_pal = iddoPal::iddo_palettes$soft_blues,
                          title = ""){
  # summarise distribution of all median preds against all covar vals
  # for each value in pfpr, find median and 95% quantile
  # might be more convenient to grab from prediction design matrix ..
  
  xax_bins <- seq(min(df[, covar]), max(df[, covar]), length.out = xax_breaks)
  df$xbinned <- cut(df[, covar], breaks = xax_bins)
  
  to_plot <- df %>% 
    group_by(xbinned) %>%
    summarise(med = median(pred),
              lower2.5 = quantile(pred, probs = c(0.025)),
              lower25 = quantile(pred, probs = c(0.25)),
              upper25 = quantile(pred, probs = c(0.75)),
              upper2.5 = quantile(pred, probs = c(0.975)),
              min = min(pred),
              max = max(pred)) %>%
    mutate(xnume = xax_bins[1:xax_breaks])
  
  p <- ggplot(data = to_plot) +
    # geom_line(aes(x = xnume, y = min)) +
    # geom_line(aes(x = xnume, y = max)) +
    geom_ribbon(aes(ymin = lower2.5, ymax = upper2.5, x = xnume, 
                    fill = "2.5% - 97.5%"),
                alpha = 0.5) +
    geom_ribbon(aes(ymin = lower25, ymax = upper25, x = xnume, 
                    fill = "25% - 75%"),
                alpha = 0.5) +
    geom_ribbon(aes(x = xnume, ymin = med, ymax = med, fill = "50%")) +
    geom_line(aes(x = xnume, y = med), col = pred_pal[1], linewidth = 1) +
    scale_fill_manual("", values = c("2.5% - 97.5%" = pred_pal[6],
                                     "25% - 75%" = pred_pal[4], 
                                     "50%" = pred_pal[1])) +
    xlab("PfPR") +
    ylab("Prevalence") +
    labs(title = title)

  if (!is.null(dat)){
    
    dat_summary <- dat %>%
      mutate(xbinned = cut(dat[, covar] %>% unlist, breaks = xax_bins)) %>%
      group_by(xbinned) %>%
      summarise(med = median(present/tested)) %>%
      mutate(xnume = xax_bins[1:nrow(.)])
    
     # p + geom_point(aes(x = pfpr, y = pred, size = tested), 
     #                data = dat, alpha = 0.6, pch = 1) +
    p + geom_point(aes(x = pfpr, y = present/tested, size = tested), data = dat,
                  pch = 1, color = "grey") +
      # killing this as it looks silly
      # geom_line(aes(x = xnume, y = med), data = dat_summary) +
      scale_size_continuous(name = "Sample size", trans = "sqrt", 
                            range = c(0.2, 5), limits = c(5, 5200), breaks = c(10, 100, 1000, 5000))
      #scale_size_continuous(trans = "sqrt", breaks = c(10, 100, 1000, 3000), "Sample size")
  }

}


library(terra)

source("code/setup.R")
source("code/build_design_matrix.R")

pfpr_unscaled <- rast("data/pfpr_rasters_afr_2025.tif")
names(pfpr_unscaled) <- paste0("pfpr_", years)

# from validation.R
mut_dat_assoc_with_preds <- lapply(names(nice_name_lookup_all), function(marker){
  extract_preds(data_path = data_path_lookup[[marker]],
                      pred_path = paste0(bb_paths[[marker]], "preds_medians.tif"),
                      buffer = BUFFER) %>%
    extract_pfpr(covs = covariates)
  
}) %>%
  setNames(names(nice_name_lookup_all)) %>%
  suppressMessages()

titlelst <- list(k13_marcse = "(a) Kelch 13", 
                 crt76 = "(b) Pfcrt K76T", 
                 mdr86 = "(c) Pfmdr1 N86Y", 
                 mdr184 = "(d) Pfmdr1 Y184F", 
                 mdr1246 = "(e) Pfmdr1 D1246Y")

plotlst <- lapply(names(nice_name_lookup_main), function(marker){
  # preds <- rast(paste0("output/", marker,"/bb_gne/preds_medians.tif"))
  # ras <- c(preds, covariates) %>%
  #   aggregate(fact = 5) # let's just make this a little more manageable
  # response_plot_runup <- response_plots_runup(ras)
  # write.csv(response_plot_runup,
  #           paste0("output/response_plot_meta/", marker, "_runup.csv"),
  #             row.names = FALSE)
  
  response_plot_runup <- read.csv(paste0("output/response_plot_meta/", marker, "_runup.csv"))
  response_plot(df = response_plot_runup,
                dat = mut_dat_assoc_with_preds[[marker]],
                xax_breaks = 100,
                title = titlelst[[marker]])
})

plotlst[[1]] + plotlst[[2]] + plotlst[[3]] + plotlst[[4]] + plotlst[[5]] +
  plot_layout(ncol = 1, guides = "collect", axis_title = "collect")
ggsave("figures/response_pfpr.png", scale = 1.5, height = 7, width = 6)



