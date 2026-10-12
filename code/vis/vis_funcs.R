# folding functions for plotting in here
library(looseVis)
library(cowplot)
library(iddoPal)
library(sf)

gg_ras_prep <- function(ras, extent = NULL, shp = NULL){
  # this was annoying me - add this to looseVis for cryin out loud
  # there has to be a better way !!
  if (!is.null(extent)){
    extent <- ext(unlist(extent))
    ras <- mask(ras, extent) %>%
      crop(extent)
    
    if (!is.null(shp)){
      # assuming we're working with sf
      shp <- st_crop(shp, extent)
    }
  }
  
  coords <- xyFromCell(ras, cells(ras))
  vals <- terra::extract(ras, coords)
  
  df <- cbind(coords, vals) %>%
    pivot_longer(starts_with("2"),
                 names_to = "lyr",
                 values_to = "val") %>%
    mutate(year = substr(lyr, 1, 4),
           tag = substr(lyr, 6, 14))
  
  list(df = df,
       shp = shp)
}



survey_effort_panel <- function(in_path, 
                                agg_factor = 1, 
                                pan = "", 
                                main = "",
                                xlab = "Longitude",
                                ylab = "Latitude",
                                lyr_names = main){
  
  test_dens <- lapply(in_path, rast) %>%
    rast()
  if (agg_factor > 1){
    test_dens <- aggregate(test_dens, agg_factor)
  }
  
  # nasty formula to convert to tests per 100kmsq (ish)
  test_dens <- test_dens * 100 / (res(test_dens)[1] * 111) ** 2
  
  test_dens <- mask(test_dens, st_as_sf(afr))
  
  surveil <- cbind(xyFromCell(test_dens, cell = cells(test_dens)), 
                   extract(test_dens, cells(test_dens)))
  names(surveil) <- c("x", "y", lyr_names)
  
  surveil <- surveil %>%
    pivot_longer(all_of(lyr_names), names_to = "lyr", values_to = "effort") %>%
    mutate(lyr = factor(lyr, levels = lyr_names))
  
  p <- ggplot() +
    geom_sf(data = st_as_sf(afr), fill = "white") + # not showing anything in the background here ...
    geom_tile(data = surveil, aes(x, y, fill = effort)) +
    geom_sf(data = st_as_sf(afr), colour = "white", fill = NA) +
    scale_fill_viridis_c(na.value = NA, bquote("Tests per ~100"~km^2), trans="sqrt") +
    facet_wrap(~ lyr) +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(paste(pan, main)) +
    theme_bw() +
    theme(legend.position = "bottom", strip.text.x = element_text(size = 12))
  
  p
}



# wrapping up plot of all-Africa posterior medians over times
pred_time_plot <- function(path, 
                           title = "",
                           pal = iddoPal::iddo_palettes$soft_blues,
                           zooms = NULL,
                           zoom_pal = NULL, 
                           alpha = 0.5,
                           show_pts = FALSE,
                           show_pt_cols = FALSE,
                           agg_fact = 10,
                           incid = NULL,
                           ylim = c(0, 1),
                           ylab = "Prevalence",
                           mask_shp = NULL,
                           cut_year = NULL){
  
  preds <- rast(paste0(path, "preds_medians.tif")) %>%
    # might be hitting a sensitivity here choosing mean ....
    # but choosing median would remove all of my peaks and not aggregating would
    # make this take 100 times longer to run :/
    aggregate(fact = agg_fact, fun = "mean", na.rm = TRUE)
  
  if (!is.null(mask_shp)){
    preds <- mask(preds, mask_shp)
  }
  
  if (!is.null(incid)){
    incid <- aggregate(incid, agg_fact, fun = "sum", na.rm = TRUE)
    keep_names <- names(preds)
    preds <- lapply(names(preds), function(x){
      yr <- str_extract(x, "^.{4}")
      idx <- which(str_extract(names(incid), ".{4}$") == yr)
      if (length(idx) == 0){idx = length(names(incid))} # assuming they're in order
      return(preds[x] * incid[[idx]])
    }) %>%
      rast()
    names(preds) <- keep_names
    ylim = c(0, max(values(preds)), na.rm=TRUE)
  }
  
  coords <- xyFromCell(preds, cells(preds))
  vals <- terra::extract(preds, coords)
  df <- cbind(coords, vals) %>%
    pivot_longer(starts_with("2"),
                 names_to = "lyr",
                 values_to = "val") %>%
    mutate(year = substr(lyr, 1, 4)) %>% # pick out year
    group_by(year) %>%
    summarise(q = list(quantile(val, c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)))) %>%
    unnest_wider(q) %>%
    ungroup() %>%
    mutate(year = as.numeric(year))
  
  if (!is.null(cut_year)){
    df <- df %>% filter(year <= cut_year)
  }
  
  message(paste0("2026 summary stats: ", 
                 df %>% filter(year == 2026) %>% dplyr::select(`50%`), "; (",
                 df %>% filter(year == 2026) %>% dplyr::select(`2.5%`), ", ",
                 df %>% filter(year == 2026) %>% dplyr::select(`97.5%`), ")"))
  
  p <- ggplot(df) +
    # removing min and max lines as they were a little distracting:
    # geom_line(aes(x = year, y = `0%`, linetype = "0% - 100%")) +
    # geom_line(aes(x = year, y = `100%`, linetype = "0% - 100%")) +
    geom_ribbon(aes(x = year, ymin = `2.5%`, ymax = `97.5%`, fill = "2.5% - 97.5%"), alpha = alpha) + #fill=pal[6]) +
    geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`, fill = "25% - 75%"), alpha = alpha) + #fill=pal[4]) +
    geom_ribbon(aes(x = year, ymin = `50%`, ymax = `50%`, fill = "50%")) + #fill=pal[1]) +
    geom_line(aes(x = year, y = `50%`), col = pal[1], linewidth = 1) +
    # scale_linetype_manual("", values = c("0% - 100%" = 2)) +
    scale_fill_manual("Quantiles", values = c("2.5% - 97.5%" = pal[6], "25% - 75%" = pal[4], "50%" = pal[1])) +
    ylab(ylab) +
    xlab("Year") +
    labs(title = title) +
    ylim(ylim[1], ylim[2]) +
    scale_x_continuous(breaks = seq(2000, 2028, 2), expand = c(0,0)) +
    geom_vline(xintercept = 2024.5, colour = iddo_palettes$BlGyRd[9], linetype = 2) +
    theme_bw() +
    theme(legend.spacing.y = unit(-10, "cm"),
          legend.background = element_rect(fill = NA))
  
  # if (!is.null(zooms)){
  #   for (i in 1:nrow(zooms)){
  #     exte <- ext(unlist(zooms[i,]))
  #     tmp <- crop(preds, exte)
  #     coords <- xyFromCell(tmp, cells(tmp))
  #     vals <- terra::extract(tmp, coords)
  #     
  #     df <- cbind(coords, vals) %>%
  #       pivot_longer(starts_with("2"),
  #                    names_to = "lyr",
  #                    values_to = "val") %>%
  #       mutate(year = substr(lyr, 1, 4),
  #              tag = substr(lyr, 11, 14)) %>% # pick out year and thingo
  #       filter(tag == "medi") %>%
  #       group_by(year) %>%
  #       summarise(med = median(val)) %>%
  #       ungroup() %>%
  #       mutate(year = as.numeric(year))
  #     
  #     p <- p + geom_line(df, mapping = aes(x = year, y = med), colour = case_pal[i])
  #   }
  #}

  if (show_pts == TRUE){
    message("Watch out! I set size limits manually!")
    mut_data <- read_rds(paste0(path, "mut_data.rds"))
    message(max(mut_data$tested))
    
    if (!is.null(mask_shp)){
      # subset to points inside mask
      mut_data <- mut_data %>%
        filter(st_as_sf(mut_data, coords = c("x", "y"), crs = st_crs(mask_shp)) %>%
                 st_intersects(st_union(mask_shp), sparse = FALSE) %>%
                 as.vector())
      
    }
    
    if (show_pt_cols == TRUE){
      mut_data <- extract_regions(mut_data)
      
      pal <- c("#c7047c", viridis(4))
      
      p <- p + geom_point(aes(x=jitter(year), y=present/tested,  
                              size=tested, col=subregion),  pch = 21,
                          mut_data) +
        scale_size_continuous(name = "Sample size", trans = "sqrt", 
                              range = c(0.2, 5), limits = c(5, 5200), breaks = c(10, 100, 1000, 5000)) +
        scale_color_manual(values = pal, "Region") +
        guides(color = guide_legend(override.aes = list(stroke = 2)))
    } else {
      p <- p + geom_point(aes(x=jitter(year), y=present/tested,  
                              size=tested), 
                          colour="grey", pch = 21,
                          mut_data) +
        scale_size_continuous(name = "Sample size", trans = "sqrt", 
                              range = c(0.2, 5), limits = c(5, 5200), breaks = c(10, 100, 1000, 5000))
    }
     # +
    # geom_boxplot(aes(x = year, y = present/tested, group = as.factor(year)),
    #              mut_data, fill = NA, outliers = FALSE)
  }
  
  p
}


summarise_ribbon <- function(path){
  preds <- rast(path) %>%
    aggregate(fact = 10)
  message("aggregating")
  
  coords <- xyFromCell(preds, cells(preds))
  vals <- terra::extract(preds, coords)
  df <- cbind(coords, vals) %>%
    pivot_longer(starts_with("2"),
                 names_to = "lyr",
                 values_to = "val") %>%
    mutate(year = substr(lyr, 1, 4)) %>% # pick out year
    group_by(year) %>%
    summarise(med = median(val)) %>%
    mutate(year = as.numeric(year))
}


extract_regions <- function(dat){
  # grab worded region for all points
  joined <- dat %>%
    st_as_sf(coords = c("x", "y"), crs = st_crs(afr)) %>%
    st_intersection(afr) %>%
    dplyr::select(name_en, subregion, tested, present, year) %>%
    suppressWarnings()
  bind_cols(st_coordinates(joined),
                   st_drop_geometry(joined))
  
}



pred_time_plot_policy <- function(path, 
                           title = "",
                           pal = iddoPal::iddo_palettes$soft_blues,
                           zooms = NULL,
                           zoom_pal = NULL, 
                           lowerupper = FALSE,
                           alpha = 0.5){
  # this one finds the median of CI bound surfaces (e.g., median of 95% lower bound surface)
  
  to_extract <- c("medians", "upper", "lower")
  to_plot <- lapply(to_extract, function(x){
    summarise_ribbon(paste0(path, "preds_", x, ".tif")) %>%
      mutate(tag = x)
  }) %>%
    do.call(what = rbind) %>%
    pivot_wider(values_from = med, names_from = tag)
    
  
  p <- ggplot(to_plot) +
    geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = "2.5% - 97.5%"), alpha = alpha) + #fill=pal[6]) +
    geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = "25% - 75%"), alpha = alpha) + #fill=pal[4]) +
    geom_ribbon(aes(x = year, ymin = medians, ymax = medians, fill = "50%")) + #fill=pal[1]) +
    geom_line(aes(x = year, y = medians), col = pal[1], linewidth = 1) +
    scale_linetype_manual("", values = c("0% - 100%" = 2)) +
    scale_fill_manual("", values = c("2.5% - 97.5%" = pal[6], "25% - 75%" = pal[4], "50%" = pal[1])) +
    ylab("Prevalence") +
    xlab("Year") +
    labs(title = title) +
    ylim(0, 1) +
    scale_x_continuous(breaks = seq(2000, 2022, 2), expand = c(0,0)) +
    theme_bw() +
    theme(legend.spacing.y = unit(-10, "cm"),
          legend.background = element_rect(fill = NA))
  p
}


pred_time_plot_crint <- function(path, 
                                  title = "",
                                  pal = iddoPal::iddo_palettes$soft_blues,
                                  zooms = NULL,
                                  zoom_pal = NULL, 
                                  alpha = 0.5,
                                  show_pts = FALSE,
                                  agg_fact = 10,
                                  ylim = NULL,
                                  proportional = FALSE,
                                  cut_year = NULL,
                                  ylab = "",
                                  plot = TRUE){
  # this could also be a fraction of annual estimated cases .........
  preds_medians <- rast(paste0(path, "preds_medians.tif"))
  yrs <- as.numeric(str_extract(names(preds_medians), "^.{4}"))
  
  preds_medians <- preds_medians %>%
    aggregate(fact = agg_fact, fun = "median", na.rm = TRUE) %>%
    sapply(function(x){median(values(x), na.rm=TRUE)})
  preds_lower <- rast(paste0(path, "preds_lower.tif")) %>%
    aggregate(fact = agg_fact, fun = "median", na.rm = TRUE) %>%
    sapply(function(x){median(values(x), na.rm=TRUE)})
  preds_upper <- rast(paste0(path, "preds_upper.tif")) %>%
    aggregate(fact = agg_fact, fun = "median", na.rm = TRUE) %>%
    sapply(function(x){median(values(x), na.rm=TRUE)})
  
  df <- data.frame(year = yrs,
                   med = preds_medians,
                   lower = preds_lower,
                   upper = preds_upper)
  
  if(!is.null(cut_year)){df <- df %>% filter(year <= cut_year)}
  
  if(is.null(ylim)){ylim <- c(min(df$lower), max(df$upper))}
  
  p <- ggplot(df) +
    geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = "2.5% - 97.5%"), alpha = alpha) +
    geom_ribbon(aes(x = year, ymin = med, ymax = med, fill = "50%")) +
    geom_line(aes(x = year, y = med), linewidth = 1, col = pal[1]) +
    scale_fill_manual("", values = c("2.5% - 97.5%" = pal[4], "50%" = pal[1])) +
    theme_bw() +
    xlab("Year") +
    labs(title = title) +
    ylab(ylab) +
    ylim(ylim[1], ylim[2]) +
    scale_x_continuous(breaks = seq(2000, 2028, 2), expand = c(0,0)) +
    geom_vline(xintercept = 2024.5, colour = iddo_palettes$BlGyRd[9], linetype = 2) +
    theme(legend.spacing.y = unit(-10, "cm"),
          legend.background = element_rect(fill = NA))
  
  p
}




factor_incidence <- function(preds, incid){
  # multiply preds by annual new cases estimate, 
  # using most recent estimate for future marker prevalence predictions
  keep_names <- names(preds)
  preds <- lapply(names(preds), function(x){
    yr <- str_extract(x, "^.{4}")
    idx <- which(str_extract(names(incid), ".{4}$") == yr)
    if (length(idx) == 0){idx = length(names(incid))} # assuming they're in order
    return(preds[x] * incid[[idx]])
  }) %>%
    rast()
  names(preds) <- keep_names
  
  preds
}


pred_time_factoring_incidence <- function(path, 
                                     incid,
                                     title = "",
                                     pal = iddoPal::iddo_palettes$soft_reds,
                                     zooms = NULL,
                                     zoom_pal = NULL, 
                                     alpha = 0.5,
                                     show_pts = FALSE,
                                     agg_fact = 10,
                                     ylim = NULL,
                                     proportional = FALSE,
                                     cut_year = NULL,
                                     ylab = "Annual cases",
                                     plot = TRUE){
  # this could also be a fraction of annual estimated cases .........
  incid <- aggregate(incid, agg_fact, fun = "sum", na.rm = TRUE)
  annual_estimated_cases <- sapply(incid, function(x){sum(values(x), na.rm = TRUE)})
  
  preds_medians <- rast(paste0(path, "preds_medians.tif"))
  yrs <- as.numeric(str_extract(names(preds_medians), "^.{4}"))
  
  preds_medians <- preds_medians %>%
    aggregate(fact = agg_fact, fun = "mean", na.rm = TRUE) %>%
    factor_incidence(incid) %>%
    sapply(function(x){sum(values(x), na.rm=TRUE)})
  preds_lower <- rast(paste0(path, "preds_lower.tif")) %>%
    aggregate(fact = agg_fact, fun = "mean", na.rm = TRUE) %>%
    factor_incidence(incid) %>%
    sapply(function(x){sum(values(x), na.rm=TRUE)})
  preds_upper <- rast(paste0(path, "preds_upper.tif")) %>%
    aggregate(fact = agg_fact, fun = "mean", na.rm = TRUE) %>%
    factor_incidence(incid) %>%
    sapply(function(x){sum(values(x), na.rm=TRUE)})
  
  df <- data.frame(year = yrs,
                   med = preds_medians,
                   lower = preds_lower,
                   upper = preds_upper)
  
  if (proportional){
    # extend estimated cases to future years - assuming we're starting in common year (2000)
    annual_estimated_cases <- c(annual_estimated_cases,
                                rep(annual_estimated_cases[length(annual_estimated_cases)],
                                    length(yrs) - length(annual_estimated_cases)))
    df <- mutate(df, across(med:upper, ~.x / annual_estimated_cases))
  }
  
  if(!is.null(cut_year)){df <- df %>% filter(year <= cut_year)}
  
  if(is.null(ylim)){ylim <- c(min(df$lower), max(df$upper))}
  
  if(plot == FALSE){
   return(paste0("2026 summary stats: ", 
                   df %>% filter(year == 2026) %>% dplyr::select(med), "; (",
                   df %>% filter(year == 2026) %>% dplyr::select(lower), ", ",
                   df %>% filter(year == 2026) %>% dplyr::select(upper), ")"))
  }
  
  message(paste0("2026 summary stats: ", 
                 df %>% filter(year == 2026) %>% dplyr::select(med), "; (",
                 df %>% filter(year == 2026) %>% dplyr::select(lower), ", ",
                 df %>% filter(year == 2026) %>% dplyr::select(upper), ")"))
  
  message(paste0("2000 summary stats: ", 
                 df %>% filter(year == 2000) %>% dplyr::select(med), "; (",
                 df %>% filter(year == 2000) %>% dplyr::select(lower), ", ",
                 df %>% filter(year == 2000) %>% dplyr::select(upper), ")"))
  
  p <- ggplot(df) +
    geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = "95% credible interval"), alpha = alpha) +
    geom_ribbon(aes(x = year, ymin = med, ymax = med, fill = "Posterior median")) +
    geom_line(aes(x = year, y = med), linewidth = 1, col = pal[1]) +
    # scale_fill_manual("", 
    #                   values = c("2.5% - 97.5%" = pal[4], "50%" = pal[1])) +
    scale_fill_manual("",
                      values = c("95% credible interval" = pal[4], 
                                 "Posterior median" = pal[1])) +
    theme_bw() +
    xlab("Year") +
    labs(title = title) +
    ylab(ylab) +
    ylim(ylim[1], ylim[2]) +
    scale_x_continuous(breaks = seq(2000, 2028, 2), expand = c(0,0)) +
    geom_vline(xintercept = 2024.5, colour = iddo_palettes$BlGyRd[1], linetype = 2) +
    theme(legend.spacing.y = unit(-10, "cm"),
          legend.background = element_rect(fill = NA))
  
  p
}




# to wrap up prediction figures for partner drug models :
map_pred_row <- function(in_path,
                         years,
                         pal,
                         field = c("50", "sd"),
                         buff = NULL,
                         exte = NULL,
                         xlab = "Longitude",
                         ylab = "Latitude",
                         legend_lim = NULL,
                         top_pan = FALSE){
  preds <- rast(in_path)
  preds <- preds[[str_extract(names(preds), "\\d{4}") %in% years]]
  
  if (!is.null(buff)){
    preds <- crop(preds, buff) %>%
      mask(buff) %>%
      trim()
  }
  
  coords <- xyFromCell(preds, cells(preds))
  vals <- terra::extract(preds, coords)
  df <- cbind(coords, vals) %>%
    pivot_longer(starts_with("2"),
                 names_to = "lyr",
                 values_to = "val") %>%
    mutate(year = substr(lyr, 1, 4)) # pick out year
  
  p <- ggplot() +
    geom_sf(data = exte, fill = NA) +
    geom_tile(data = df, 
              mapping = aes(x = x, y = y, fill = val)) +
    geom_sf(data = exte, fill = NA, col = "grey") +
    facet_wrap(~year, nrow = 1) +
    #xlab(xlab) +
    ylab(ylab) +
    #scale_x_continuous(breaks = seq(-20, 40, 20)) +
    #scale_y_continuous(breaks = seq(-20, 40, 20)) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_blank(),
          legend.justification = "left",
          axis.title.x = element_blank(),
          #axis.title.y = element_text(angle = 0, hjust = 1),
          #axis.title.y = element_blank(),
          #title = element_blank(),
          panel.spacing = unit(0, "lines"))
  
  if (!top_pan){
    p <- p + theme(strip.background = element_blank(),
                   strip.text.x = element_blank())
  }
  
  if (field == "50"){
    if (!is.null(legend_lim)){
      p <- p + scale_fill_gradientn(name = "Prevalence",
                                    colors = pal, 
                                    breaks = c(0, 0.5, 1), 
                                    labels = c("0  (all wildtype)", "0.5", "1  (all mutant)"),
                                    limits = legend_lim)
    } else {
      p <- p + scale_fill_gradientn(name = "Prevalence",
                                    colors = pal, 
                                    limits = legend_lim)
    }
  } else {
    p <- p + scale_fill_gradientn(colors = pal, 
                                  name = "Estimate SD",
                                  limits = legend_lim)
  }
  
  p
}


# convert to gg
# add option to facet into time windows
# obs_prev_panel_base <- function(data_path, 
#                                 pred_path, 
#                                 main = "", 
#                                 show_nas = FALSE, 
#                                 #pal = colorRamp(viridis(10)), 
#                                 pal = colorRamp(iddo_palettes$BlGyRd),
#                                 xlim = c(0,1), ylim = c(0,1)){
#   mut_data <- setup_mut_data(data_path, min_year = MIN_YEAR)
#   preds <- rast(pred_path)
#   yrs_pred <- str_extract(names(preds), "\\d{4}")
#   
#   mut_data$pred <- NA
#   yrs_to_extract <- unique(mut_data$year)
#   for (yr in yrs_to_extract){
#     if (yr %in% yrs_pred){
#       idx <- which(mut_data$year == yr)
#       val <- terra::extract(preds[[paste0(yr, "_post_median")]], 
#                             mut_data[idx, c("x", "y")],
#                             ID = FALSE)
#       mut_data[idx, "pred"] <- val
#     }
#   }
#   
#   mut_data = mut_data[!is.na(mut_data$pred),]
#   # mut_data$diffs = abs(mut_data$present / mut_data$tested - mut_data$pred)
#   mut_data$diffs = mut_data$present / mut_data$tested - mut_data$pred
#   message(paste(min(mut_data$diffs), max(mut_data$diffs)))
#   mut_data <- arrange(mut_data, abs(diffs))
#   mut_data$diffs = mut_data$diffs / 2 + 0.5 # hopefully grey ends up where diffs == 0?
#   message(paste(min(mut_data$diffs), max(mut_data$diffs)))
#   
#   plot(mut_data$present / mut_data$tested, 
#        mut_data$pred, cex = cex_transform(mut_data$tested) * 4,
#        xlab = "Observed prevalence", ylab = "Predicted prevalence", 
#        xlim = xlim, ylim = ylim,
#        col = rgb(pal(mut_data$diffs), maxColorValue = 255))
#   abline(a = 0, b = 1)
#   
#   plot(st_geometry(afr))
#   points(mut_data$x, mut_data$y, pch = 16,
#          col = rgb(pal(mut_data$diffs), maxColorValue = 255))
#   
#   mtext(outer = TRUE, text = main)
# }
# here's an e.g.:
# obs_prev_panel_base("data/clean/moldm_k13_nomarker.csv",
#                "output/k13/gneiting_sparse/preds_all.grd",
#                "k13 gneiting", xlim = c(0, 0.4), ylim = c(0, 0.4))



coverages_inner <- function(path){
  # wrapping extraction for a single model
  coverages <- read.csv(paste0(path, "coverages/all_coverages.csv"))
  covs_all <- coverages[seq(1, nrow(coverages), 2),] %>%
    colSums(na.rm=TRUE)
  covs_pos <- coverages[seq(2, nrow(coverages), 2),] %>%
    colSums(na.rm=TRUE)
  
  covs_all <- covs_all[grep("X", names(covs_all))] / covs_all["n_landed"]
  covs_pos <- covs_pos[grep("X", names(covs_pos))] / covs_pos["n_non_zero"]
  widths <- gsub("X", "", names(covs_all))[grep("X", names(covs_all))] %>% as.numeric()
  
  dat <- data.frame(w = widths / 100,
                    pa = covs_all,
                    po = covs_pos) %>%
    pivot_longer(cols = c(pa, po), names_to = "group", values_to = "value") %>%
    mutate(group = ifelse(group == "pa", "All points", "Presences only"))
  
  if (!grepl("k13", path)){
    dat <- dat %>% filter(group != "Presences only")
    # dat <- dat %>% dplyr::select(-c(group))
  }
  dat
}


coverages_fig <- function(path){
  # wish I had included years here ... wait a sec I did
  # yrs <- concat_coverages(path)
  
  # removing earlier code when I was comparing binom and betabinom models:
  # if (length(path) > 1){
  #   dat <- lapply(path, function(p){
  #     coverages_inner(p) %>% 
  #       mutate(mod = str_extract(p, "(?<=/)[^/]+(?=/[^/]*$)"))
  #     }) %>%
  #     do.call(what = rbind) %>%
  #     mutate(mod = case_when(mod == "bb_gne" ~ "Beta-binomial",
  #                            mod == "gneiting_sparse" ~ "Binomial",
  #                            TRUE ~ mod))
  # } else {
  #   dat <- coverages_inner(path)
  # }
  # 
  # # plot(widths, covs_pos, type="l", 
  # #      ylim = c(0, 100), xlim = c(0, 100))
  # # lines(widths, covs_all, type = "l")
  # # abline(a = 0, b = 1)
  # 
  # if ("mod" %in% names(dat)){
  #   p <- ggplot(dat) +
  #     geom_line(aes(x = w, y = value, color = mod, linetype = group)) +
  #     geom_abline(slope = 1, intercept = 0, col="darkgrey") +
  #     scale_color_manual(values = iddo_palettes$iddo, name="Model") +
  #     scale_linetype("")
  # } else {
  #   p <- ggplot(dat) +
  #     geom_line(aes(x = w, y = value, color = group)) +
  #     geom_abline(slope = 1, intercept = 0, col="darkgrey") +
  #     scale_color_manual(values = iddo_palettes$iddo, name="")
  # }
  
  dat <- lapply(path, function(p){
    message(p)
    message(str_extract(p, "(?<=/)[^/]+(?=/)"))
    coverages_inner(p) %>%
      mutate(marker = str_extract(p, "(?<=/)[^/]+(?=/)"),
             marker = factor(unlist(nice_name_lookup_all[marker]),
                                levels = nice_name_lookup_all))
  }) %>%
    do.call(what = rbind)
  
  p <- ggplot(dat) +
    geom_line(aes(x = w, y = value, color = marker, linetype = group)) +
    geom_abline(slope = 1, intercept = 0, col = "darkgrey") +
    scale_color_manual(values = viridis(5), name = "Marker")
  
  p + 
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    xlab("Credible interval width") +
    ylab("Coverage probability") +
    theme_bw()
}




# unfinished: ...
# plot_partner_markers <- function(shp, buff, 
#                                  markers = c("mdr86", "mdr184", "mdr1246", "crt76")){
#   for (marker in markers){
#     mut_data <- setup_mut_data(get_input_dir(marker), min_year = MIN_YEAR)
#   }
#   
# }

thresholds_fig <- function(path,
                           years_to_plot,
                           thresholds,
                           out = ""){
  preds <- rast(paste0(path, "preds_medians.tif"))
  preds <- preds[[str_extract(names(preds), "\\d{4}") %in% years_to_plot]] #%>%
  #  aggregate(fact = 10)
  
  df <- gg_ras_prep(preds)$df

  tmp <- df %>%
    mutate(val = cut(val, breaks = c(0, thresholds, 100)))
  
  pal <- colorRampPalette(iddo_palettes$Blues)(length(thresholds) + 1)
  
  ggplot(tmp) +
    geom_sf(data = afr, col = "grey70") +
    geom_tile(aes(x = x, y = y, fill = val)) +
    scale_fill_manual(values = rev(pal), "Median Kelch 13\nprevalence") +
    facet_wrap(~year) +
    geom_sf(data = afr, col = "grey40", fill = NA) +
    theme_bw() +
    labs(x = "Longitude", y = "Latitude") +
    theme(legend.position = "bottom")
  
  if (out != ""){
    ggsave(paste0("figures/thresholds_", out, ".png"),
           height = 5, width = 10)
  }
}

# path <- "output/k13_marcse/gneiting_sparse/"
# years_to_plot <- c("2012", "2018", "2024")
# thresholds <- c(0.02, 0.05, 0.1, 0.25)
# thresholds_fig("output/k13_marcse/gneiting_sparse/",
#                c("2012", "2018", "2024"),
#                c(0.02, 0.05, 0.1, 0.25),
#                "k13m_bi")
# thresholds_fig("output/k13_marcse/bb_gne/",
#                c("2012", "2018", "2024"),
#                c(0.02, 0.05, 0.1, 0.25),
#                "k13m_bb")


lower_upper_panel_crints <- function(path, 
                                     main = "", 
                                     show_nas = FALSE, 
                                     pal = colorRamp(iddo_palettes$BlGyRd),
                                     xlim = c(0,1), # define limits to pred/obs panel
                                     ylim = c(0,1), # define limits to pred/obs panel
                                     facet_bins = NULL, # apply facets over time?
                                     ave_tag = "_50", # mean? median? what are the surfaces called in the stack?
                                     buffer = 1, # option to reland points?
                                     bb = NULL){
  # assess credible interval coverage of observed marker prevalences
  mut_data <- read_rds(paste0(path, "mut_data.rds")) %>%
    arrange(present/tested)
  lower <- rast(paste0(path, "preds_lower.tif"))
  upper <- rast(paste0(path, "preds_upper.tif"))
  yrs_pred <- str_extract(names(upper), "\\d{4}")
  
  # get predictions for each row in `mut_data`
  mut_data$lower <- NA
  mut_data$upper <- NA
  yrs_to_extract <- unique(mut_data$year)
  for (yr in yrs_to_extract){
    if (yr %in% yrs_pred){
      idx <- which(mut_data$year == yr)
      vall <- terra::extract(lower[[paste0(yr, "_2.5")]], 
                             mut_data[idx, c("x", "y")],
                             ID = FALSE, search_radius = buffer)
      valu <- terra::extract(upper[[paste0(yr, "_97.5")]], 
                             mut_data[idx, c("x", "y")],
                             ID = FALSE, search_radius = buffer)
      if(ncol(vall) < 3){
        # idk why we have to have an inconsistent return when |idx| == 1
        mut_data[idx, "lower"] <- vall[1,1]
        mut_data[idx, "upper"] <- valu[1,1]
      } else{
        mut_data[idx, "lower"] <- vall[, paste0(yr, "_2.5")]
        mut_data[idx, "upper"] <- valu[, paste0(yr, "_97.5")]
      }
    }
  }
  
  message(nrow(mut_data))
  mut_data <- mut_data %>% 
    filter(!is.na(lower) & !is.na(upper)) %>%
    mutate(idx = 1:nrow(.), 
           covered = case_when(present/tested < lower ~ "Over-prediction",
                               present/tested > upper ~ "Under-prediction",
                               present/tested > lower & present/tested < upper ~ "Covered",
                               TRUE ~ NA))
  #covered = ifelse(present/tested > lower & present/tested < upper, TRUE, FALSE))
  message(nrow(mut_data))
  
  mut_data_long <- mut_data %>%
    pivot_longer(cols = c(lower, upper), names_to = "lu", values_to = "value")
  
  ggplot(mut_data) +
    geom_point(aes(x = idx, y = present/tested,
                   size = tested, 
                   col = covered), 
               pch = 1) +
    geom_line(data = mut_data_long, aes(x = idx, y = value, group = idx, col = covered)) +
    scale_colour_manual(values = iddo_palettes$BlGyRd[c(5,2,8)], "CI covers\nobservation") +
    scale_size_continuous(range = c(0.1, 6), "Tested") +
    xlab("Index") +
    ylab("Observed prevalence") +
    scale_x_continuous(expand = c(0.01,0.01)) + # remove excess padding
    theme_bw()
  # for some reason, this legend.position doesn't play with get_legend
  # theme(legend.position = "bottom")
  
}


