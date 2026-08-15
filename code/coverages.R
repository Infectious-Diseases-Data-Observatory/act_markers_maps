suppressWarnings(suppressMessages(source("code/setup.R")))
suppressMessages(source("code/build_design_matrix.R"))
source("code/predict_to_raster.R")
source("code/betabinomial_p_rho.R")
source("code/validation_funcs.R") # helpers are mostly in here

set.seed(0748)

# first, take posterior predictions
# sapply(names(nice_name_lookup_all), function(marker){
#   message(marker)
#   post_pred_sims_please(marker)
# })

# second, pass posterior predictions (prevalence) through obs model (counts)
# and compare quantiles of samples to observed number of cases with marker
sim_coverages <- lapply(names(nice_name_lookup_all), 
                        function(marker){
  message(marker)
  coverage_probabilities_from_observation_model(marker,
                                                probs = seq(0, 1, 0.01),
                                                nsim = 20, # this is just for rho ..
                                                new_samps = FALSE) # can turn this off now I have samps
})
names(sim_coverages) <- names(nice_name_lookup_all)

# have a geez at positives only for kelch
sim_coverages_pos_only <- lapply(names(nice_name_lookup_all)[grepl("k13", names(nice_name_lookup_all))], 
     function(marker){
       message(marker)
       coverage_probabilities_from_observation_model(marker,
                                                     pos_only = TRUE,
                                                     probs = seq(0, 1, 0.01),
                                                     nsim = 20, # this is just for rho ..
                                                     new_samps = FALSE)
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

write.csv(sim_coverages, "output/sim_coverages.csv", row.names = FALSE)

library(viridisLite)
library(cowplot)
kelch_pal <- c("#14B1E7", viridis(5))
partner_pal <- c("#c7047c", "#7ECE7E", "#174D97", "#E37210")


p1 <- ggplot(data = sim_coverages %>%
               mutate(markerf = factor(unlist(nice_name_lookup_all[marker]),
                                       levels = nice_name_lookup_all)) %>%
               filter(grepl("k13", marker)),
             aes(x = widths, y = cover)) +
  geom_line(aes(group = interaction(markerf, recs), col = markerf, linetype = recs)) +
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
  geom_line(aes(group = interaction(markerf, recs), col = markerf, linetype = recs)) +
  geom_abline(slope = 1, col = "grey") +
  scale_color_manual(values = partner_pal) + #, drop = FALSE) +
  xlab("Posterior predictive interval width") +
  ylab("Coverage probability") +
  scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  theme_bw()

p <- plot_grid(p1, p2)
ggsave("figures/coverages.png", p, height = 5, width = 9)


# now for PIT-ECDFs
pit_ecdfs_out <- lapply(names(nice_name_lookup_all), function(marker){
  message(marker)
  pit_ecdfs(marker, new_samps = FALSE)
})
names(pit_ecdfs_out) <- names(nice_name_lookup_all)

pit_ecdfs_out <- lapply(names(pit_ecdfs_out), function(x){
  mutate(pit_ecdfs_out[[x]], marker = x)
}) %>%
  do.call(what = rbind)



p3 <- ggplot(bind_rows(pit_ecdfs_out %>% 
                         filter(name == "leq" & 
                                  grepl("k13", marker)) %>%
                         mutate(dat = "All records"), 
                       pit_ecdfs_out %>% 
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

p4 <- ggplot(pit_ecdfs_out %>% 
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

p <- plot_grid(p3, p4)
ggsave("figures/pit_ecdfs.png", p, height = 5, width = 9)

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
ggsave("figures/coverages_all.png", width = 9, height = 9)




