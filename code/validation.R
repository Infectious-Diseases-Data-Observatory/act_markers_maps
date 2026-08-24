# script for validation outputs

# # for nearest neighbour index - from Foo and Flegg
# library(tensorflow)
source("code/setup.R")
source("code/build_design_matrix.R") # for year scaling

# bringing in some backend functions from greta.gp:
# source("~/greta.gp.st.on.earth/R/tf_kernels.R")
source("code/betabinomial_p_rho.R")

library(looseVis)
library(iddoPal)
library(cowplot)

source("code/validation_funcs.R")

val_table_row_order <- data.frame(marker = names(nice_name_lookup_all),
                                  ord = c(1, 7:10, 2:6))


# read in mut_data and associate each record with predicted prevalence for 
# relevant year:
mut_dat_assoc_with_preds <- lapply(names(nice_name_lookup_all), function(marker){
  extract_preds(data_path = data_path_lookup[[marker]],
                pred_path = paste0(bb_paths[[marker]], "preds_medians.tif"),
                buffer = BUFFER) # 5000 is sufficient to drag all points onto mask .. although what pfpr did I assign them during fitting?
}) %>%
  setNames(names(nice_name_lookup_all)) %>%
  suppressMessages()


# now we're ready to work with cv preds - here are the preds that I've slurmed for 
# all points
mut_dat_assoc_with_preds_cv <- lapply(names(nice_name_lookup_all), function(marker){
  read.csv(paste0("output/", marker, "/bb_gne/mut_dat_cv_preds_extracted.csv"))
}) %>%
  setNames(names(nice_name_lookup_all))

cv_folds <- lapply(names(nice_name_lookup_all), function(marker){
  read_rds(paste0("output/", marker, "/bb_gne/cv_folds_spat.rds"))
}) %>%
  setNames(names(nice_name_lookup_all))


# grab LFO-CV preds also - main models only as SNP models were tricky
mut_dat_assoc_with_preds_lfo <- lapply(names(nice_name_lookup_main),
                                       function(marker){
   read.csv(paste0("output/", marker, "/bb_gne/lfo/mut_dat_cv_preds_extracted.csv"))
}) %>%
  setNames(names(nice_name_lookup_main))

# this will be in reference to a different number of records compared to 
#what's in mut_dat_assoc_with_preds_lfo
# lfo_folds <- lapply(names(nice_name_lookup_main), function(marker){
#   read_rds(paste0("output/", marker, "/bb_gne/lfo/cv_folds_lfo.rds"))
# })


# baseline predictions: just assign annual marker median
mut_dat_preds_baseline <- lapply(names(nice_name_lookup_all), function(marker){
  baseline_preds(mut_dat_assoc_with_preds[[marker]])
}) %>%
  setNames(names(nice_name_lookup_all))
# might be prudent to trim regions I'm looking at for SNP models ...need to think further about that


# summarise
rmses <- lapply(mut_dat_assoc_with_preds, function(x){rmse(x)})
rsq <- lapply(mut_dat_assoc_with_preds, function(x){unadjusted_rsq(x)})

# crack into rsqs for held-out models and add n folds included in mean
cv_stats <- lapply(names(nice_name_lookup_all), 
                   function(x){cv_val(mut_dat_assoc_with_preds_cv[[x]])}) %>%
  setNames(names(nice_name_lookup_all))

cv_stats_lfo <- lapply(names(nice_name_lookup_main), 
                   function(x){cv_val(mut_dat_assoc_with_preds_lfo[[x]])}) %>%
  setNames(names(nice_name_lookup_main))

rmse_baseline <- lapply(mut_dat_preds_baseline, function(x){rmse(x)})
rsq_baseline <- lapply(mut_dat_preds_baseline, function(x){unadjusted_rsq(x)})

cv_summary <- lapply(cv_stats, function(x){
  x[-c(1)] %>% # remove foldwise stats
    as.data.frame()
}) %>%
  do.call(what = rbind) %>%
  mutate(marker = rownames(.))

cv_summary_lfo <- lapply(cv_stats_lfo, function(x){
  x[-c(1)] %>% # remove foldwise stats
    as.data.frame()
}) %>%
  do.call(what = rbind) %>%
  mutate(marker = rownames(.)) %>%
  rename_with(.cols = -c(marker), ~paste0(.x, "_lfo"))

library(xtable)
dat <- data.frame(mod = unlist(nice_name_lookup_all[names(rmses)]),
                  base_rmse = unlist(rmse_baseline),
                  base_rsq = unlist(rsq_baseline),
                  rmse = unlist(rmses),
                  rsq = unlist(rsq)) %>%
  mutate(marker = rownames(.)) %>%
  left_join(cv_summary) %>%
  left_join(cv_summary_lfo) %>%
  left_join(val_table_row_order) %>%
  mutate(across(matches("rsq|rmse"), ~ sprintf("%.3f", .x))) %>%
  mutate(across(ends_with("lfo"), ~ ifelse(.x == "NA", "", .x))) %>%
  mutate(across(starts_with("base"), ~ ifelse(.x == "NA", "*", .x))) %>%
  arrange(ord) %>%
  dplyr::select(-c(marker, ord))

dat

colnames(dat) <- c("", "RMSE", "$r^2$", "RMSE", "$r^2$", "Mean", "SD", "Mean", 
                   "SD", "n", "Mean", "SD", "Mean", "SD", "n")
tab <- xtable(dat)
align(tab) <- c(rep("c", 2), "|", rep("c", 2), "|", rep("c", 2), "|", rep("c", 5),
                "|", rep("c", 5))
addtorow <- list()
addtorow$pos <- list(-1, -1, -1)
addtorow$command <- c("\\hline & \\multicolumn{2}{c|}{Baseline} & \\multicolumn{11}{c}{Spatiotemporal GP} \\\\\n \\cline{4-15}",
                      " & & & \\multicolumn{2}{c|}{Full dataset} & \\multicolumn{5}{c|}{10-fold SS-CV} & \\multicolumn{5}{c}{LFO-CV (2019--2024)} \\\\\n", # \\cline{6-9}
                      "& & & & & \\multicolumn{2}{c}{RMSE} & \\multicolumn{2}{c}{$r^2$} & & \\multicolumn{2}{c}{RMSE} & \\multicolumn{2}{c}{$r^2$} & \\\\\n")
print(tab, 
      sanitize.text.function=function(x){x}, 
      include.rownames = FALSE,
      #include.colnames = FALSE,
      hline.after = c(0, 10),
      add.to.row = addtorow)

################################################################################
# check distribution of observed prevalences and how that varies through blocks ..


pal <- c("#14B1E7", "#440154FF", "#135ced", "#5DC863FF", "#FDE725FF", "#c7047c", 
  "#E37210","#f2bbee", "#9d13ed", "#13edde")

ggplot(data = mut_dat_assoc_with_preds_cv$k13_marcse %>%
         mutate(fold = as.factor(fold))) +
  geom_sf(data = afr) +
  geom_point(aes(x = x, y = y, col = fold)) +
  scale_color_discrete("Fold") +
  theme(axis.title = element_blank())
ggsave("figures/blocks_eg.png", height = 7, width = 7)

ggplot(data = mut_dat_assoc_with_preds_cv$k13snp_A675V %>%
         mutate(fold = as.factor(fold))) +
  geom_histogram(aes(fill = fold, x = present/tested)) +
  facet_wrap(~fold, ncol = 1, scales = "free_y") +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) +
  scale_y_continuous(trans = "log")

# I could do this for all of the models ... but that might be too much

################################################################################
# plots of residuals

abcde = c("a", "b", "c", "d", "e")
tmp = lapply(mut_dat_assoc_with_preds[grepl("k13", names(mut_dat_assoc_with_preds))], 
             obs_prev_panel, legend_position = "bottom")

# splitting the above in three:
p <- plot_grid(tmp$k13_marcse,
               tmp$k13snp_A675V,
               tmp$k13snp_C469Y,
               ncol = 1) +
  theme(plot.margin = margin(0.7, 0, 0, 0, unit = "cm"))
ggsave("figures/obs_prev_kelch_1.png", 
       p + geom_text(aes(x = 0, 
                         y = rev(seq(0.34, 1.01, length.out = 3)), 
                         label = paste0("(", abcde[1:3], ") ", nice_name_lookup_all[c(1, 6, 7)])),
                     hjust = 0), 
       height = 7.5, scale = 1.5, width = 6)

p <- plot_grid(tmp$k13snp_P441L,
               tmp$k13snp_R561H,
               tmp$k13snp_R622I,
               ncol = 1) +
  theme(plot.margin = margin(0.7, 0, 0, 0, unit = "cm"))
ggsave("figures/obs_prev_kelch_2.png", 
       p + geom_text(aes(x = 0, 
                         y = rev(seq(0.34, 1.01, length.out = 3)), 
                         label = paste0("(", abcde[1:3], ") ", nice_name_lookup_all[8:10])),
                     hjust = 0), 
       height = 7.5, scale = 1.5, width = 6)

tmp = lapply(mut_dat_assoc_with_preds[!grepl("k13", names(mut_dat_assoc_with_preds))], 
             obs_prev_panel, legend_position = "left")

p <- plot_grid(tmp$crt76, 
               tmp$mdr86,
               tmp$mdr184,
               tmp$mdr1246, 
               ncol = 1) +
  theme(plot.margin = margin(0.2, 0, 0, 0, unit = "cm"))
ggsave("figures/obs_prev_partners.png", 
       p + geom_text(aes(x = 0, 
                         y = rev(seq(0.245, 1, length.out = 4)), 
                         label = paste0("(", abcde[1:4], ") ", nice_name_lookup_all[2:5])),
                     hjust = 0), 
       height = 7.5, scale = 1.5, width = 6)


# tmp2 <- mut_dat_assoc_with_preds$mdr1246 %>%
#   mutate(diff = present/tested - pred) %>%
#   # filter(diff > 0.5)
#   filter(y > -5 & y < 5 & x > 28 & y < 37)
# 
# # it looks as though we have lots of points on top of each other .......
# ggplot() +
#   geom_sf(data = afr %>% filter(name %in% c("Uganda", "Kenya"))) +
#   geom_jitter(aes(x = x, y = y, size = tested, col = diff),
#               height = 0.1,
#               width = 0.1,
#               alpha = 0.5,
#              tmp2 %>%
#                mutate(year_bin = cut(year, c(1999, 2010, 2019, 2020, 2025))) %>%
#                arrange(diff) %>%
#                mutate(abs_diff = abs(diff))) +
#   facet_wrap(~year_bin) +
#   scale_color_viridis_c(option = "H")


# ggplot() +
#   geom_sf(data = afr %>% filter(name %in% c("Uganda", "Kenya"))) +
#   geom_jitter(aes(x = x, y = y, size = tested, col = abs_diff),
#               height = 0.1,
#               width = 0.1,
#               alpha = 0.5,
#               tmp2 %>%
#                 mutate(year_bin = cut(year, c(1999, 2010, 2015, 2020, 2025))) %>%
#                 arrange(diff) %>%
#                 mutate(abs_diff = abs(diff))) +
#   facet_wrap(~year_bin) +
#   scale_color_viridis_c()



# given predicted prevalence, take posterior samples at location of all observations
# and compare quantiles of samples to observed number of cases with marker
sim_coverages <- lapply(names(nice_name_lookup_all), function(marker){
  message(marker)
  coverage_probabilities_from_observation_model(mut_dat_assoc_with_preds[[marker]],
                                                bb_paths[[marker]],
                                                probs = seq(0, 1, 0.01),
                                                nsim = 500)
})
# (all of those samples and their summarisation takes a bit of a while)
names(sim_coverages) <- names(nice_name_lookup_all)

sim_coverages_pos_only <- lapply(names(nice_name_lookup_all)[grepl("k13", names(nice_name_lookup_all))], 
                                 function(marker){
  message(marker)
  coverage_probabilities_from_observation_model(mut_dat_assoc_with_preds[[marker]] %>%
                                                  filter(present > 0),
                                                bb_paths[[marker]],
                                                probs = seq(0, 1, 0.01),
                                                nsim = 500)
})
names(sim_coverages_pos_only) <- names(nice_name_lookup_all)[grepl("k13", names(nice_name_lookup_all))]

sim_coverages <- bind_rows(
  lapply(names(sim_coverages), function(x){
    mutate(sim_coverages[[x]], marker = x)
  }) %>%
  do.call(what = rbind) %>%
  mutate(recs = "All records"),
  lapply(names(sim_coverages_pos_only), function(x){
    mutate(sim_coverages_pos_only[[x]], marker = x)
  }) %>%
    do.call(what = rbind) %>%
    mutate(recs = "Presences only")
)




posterior_predictive_ecdfs <- lapply(names(nice_name_lookup_all), function(marker){
  message(marker)
  posterior_predictive_check(mut_dat_assoc_with_preds[[marker]],
                             bb_paths[[marker]])
})
names(posterior_predictive_ecdfs) <- names(nice_name_lookup_all)

posterior_predictive_ecdfs <- lapply(names(posterior_predictive_ecdfs), function(x){
  mutate(posterior_predictive_ecdfs[[x]], marker = x)
}) %>%
  do.call(what = rbind)

# have ended up excluding:
# p1 <- coverages_fig(paste0("output/", c("k13_marcse", "crt76", "mdr86", 
#                                         "mdr184", "mdr1246"), "/bb_gne/"))

iddo_palettes$iddo
pal <- c(  viridis(6), "#E37210")
kelch_pal <- c(iddoblue, viridis(5))
partner_pal <- c("#c7047c", iddo_palettes$iddo_new[c(1, 3, 5)])
partner_pal <- c("#c7047c", "#7ECE7E", "#174D97", "#E37210")


p1 <- ggplot(data = sim_coverages %>%
               mutate(markerf = factor(unlist(nice_name_lookup_all[marker]),
                                      levels = nice_name_lookup_all)) %>%
               filter(grepl("k13", marker)),
             aes(x = widths, y = cover)) +
  # geom_point(aes(group = interaction(marker, recs), col = marker), size = 0.5) +
  geom_line(aes(group = interaction(markerf, recs), col = markerf, linetype = recs)) + #,
           # show.legend = TRUE) +
  geom_abline(slope = 1, col = "grey") +
  scale_color_manual(values = kelch_pal) + #, drop = FALSE) +
  scale_linetype_manual(values = c("solid", "longdash")) +
  xlab("Posterior predictive interval width") +
  ylab("Coverage probability") +
  scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  theme_bw()

p2 <- ggplot(data = sim_coverages %>%
               mutate(markerf = factor(unlist(nice_name_lookup_all[marker]),
                                      levels = nice_name_lookup_all)) %>%
               filter(!grepl("k13", marker)), 
             aes(x = widths, y = cover)) +
  # geom_point(aes(group = interaction(marker, recs), col = marker), size = 0.5) +
  geom_line(aes(group = interaction(markerf, recs), col = markerf, linetype = recs)) + #,
            # show.legend = TRUE) +
  geom_abline(slope = 1, col = "grey") +
  scale_color_manual(values = partner_pal) + #, drop = FALSE) +
  xlab("Posterior predictive interval width") +
  ylab("Coverage probability") +
  scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  theme_bw()

p3 <- ggplot(bind_rows(posterior_predictive_ecdfs %>% 
                         filter(name == "leq" & 
                                 grepl("k13", marker)) %>%
                         mutate(dat = "All records"), 
                       posterior_predictive_ecdfs %>% 
                         filter(name == "leq" & 
                                  grepl("k13", marker) &
                                  #marker == "k13_marcse" & 
                                  present != 0) %>%
                         mutate(dat = "Presences only")) %>%
               mutate(marker = factor(unlist(nice_name_lookup_all[marker]),
                                      levels = nice_name_lookup_all)),
             aes(x = value)) +
  stat_ecdf(geom = "step", 
            aes(group = interaction(marker, dat), 
                col = marker, linetype = dat)) + #, show.legend = TRUE) +
  geom_abline(intercept = 0, slope = 1, col = "grey") +
  scale_color_manual("Kelch 13 markers", values = kelch_pal) + #, drop = FALSE) + 
  scale_linetype_manual("", values = c("solid", "longdash")) +
  ylab("ECDF") +
  xlab("Pr(posterior predictive samples <= observed data)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(legend.box = "horizontal")

p4 <- ggplot(posterior_predictive_ecdfs %>% 
                         filter(name == "leq" & !grepl("k13", marker)) %>%
                         mutate(dat = "All records") %>%
               mutate(marker = factor(unlist(nice_name_lookup_all[marker]),
                                      levels = nice_name_lookup_all)),
             aes(x = value)) +
  stat_ecdf(geom = "step", 
            aes(group = interaction(marker, dat), 
                col = marker, linetype = dat)) +
  geom_abline(intercept = 0, slope = 1, col = "grey") +
  scale_color_manual("Partner drug markers", values = partner_pal) + 
  scale_linetype("", guide = "none") +
  ylab("ECDF") +
  xlab("Pr(posterior predictive samples <= observed data)") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw()
# Kelch 13: high number of zeroes and low probabilities: 100% of samples leq observation

# leg = get_legend(p2)
# plot_grid(
#   plot_grid(# p1 + theme(legend.position = "none"),
#             p1 + theme(legend.position = "none"), 
#             p2 + theme(legend.position = "none"), 
#             p3 + theme(legend.position = "none"), 
#             p4 + theme(legend.position = "none"), 
#             ncol = 1, align = "v",
#             labels = c("(a)", "(b)"), label_fontface = "plain", label_x = 0.11, label_y = 0.98),
#   leg, nrow = 1, rel_widths = c(1, 0.35)
# )
# # ggsave("figures/coverages_with_unsim.png", width = 7, height = 7)
# ggsave("figures/coverages.png", width = 9, height = 9)

pan_mar <- c(0.2,0.3,0.2,0.3)
leg1 = get_legend(p3 + theme(legend.justification = c(0.5,1)))
leg2 = get_legend(p4 + theme(legend.justification = c(0.5,1)))
plot_grid(
  plot_grid(# p1 + theme(legend.position = "none"),
    p1 + theme(legend.position = "none",
               plot.margin = unit(pan_mar, "cm")), 
    p2 + theme(legend.position = "none",
               plot.margin = unit(pan_mar, "cm")), 
    p3 + theme(legend.position = "none",
               plot.margin = unit(pan_mar, "cm")), 
    p4 + theme(legend.position = "none",
               plot.margin = unit(pan_mar, "cm")), 
    ncol = 2, align = "v",
    labels = c("(a)", "(b)", "(c)", "(d)"), label_fontface = "plain", 
    label_x = -0.01, label_y = 0.96),
  plot_grid(leg1, leg2, ncol = 2, align = "h", axis = "t"), 
  nrow = 2, rel_heights = c(1, 0.27)
)
ggsave("figures/coverages.png", width = 9, height = 9)


# mut_data <- mut_data %>%
#   mutate(nn = nn_measure(mut_data, 
#                          draws_path),
#          nnplot = 2**nn,
#          nnplot = nnplot/max(nnplot))
# 
# ggplot(mut_data) +
#   geom_point(aes(x = x, y = y, col = nnplot))
# 
# ggplot(mut_data) +
#   geom_point(aes(x = present/tested, y = abs(pred - present/tested), col = nnplot)) +
#   scale_color_viridis_c("Mean distance\nto other points")
# 
# ggplot(mut_data) +
#   geom_sf(data = afr) +
#   geom_point(aes(x = x, y = y, col = nnplot), alpha = 0.3) +
#   scale_color_viridis_c("Mean distance\nto other points") +
#   xlab("Longitude") +
#   ylab("Latitude")

# # e.g.:
# # might want to re-land some points inside of model fitting
obs_prev_panel("data/clean/moldm_marcse_k13_nomarker.csv",
               "output/k13_marcse/bb_gne/preds_medians.tif",
               xlim = c(0, 0.6), ylim = c(0, 0.6),
               ave_tag = "_50", buffer = 100000, bb = c(27, 37, -5,  5))

tmp <- lapply(names(data_path_lookup), function(marker){
  message(marker)
  obs_prev_panel(data_path_lookup[[marker]],
                 paste0(bb_paths[[marker]], "preds_medians.tif"),
                 xlim = c(0, 1), ylim = c(0, 1),
                 ave_tag = "_50", buffer = 100000, bb = c(27, 37, -5,  5))
})
names(tmp) <- names(data_path_lookup)

# grid these into supp:
p <- plot_grid(tmp$k13_marcse, tmp$crt76, tmp$mdr86, tmp$mdr184, tmp$mdr1246,
          nrow = 3, ncol = 2)
ggsave("figures/obs_prev_panelled.png", p, height = 9, width = 8, scale = 1.4)
# legends will need a fiddle
# could be 3 * 5 with UGA inset ...





