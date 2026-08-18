# library(bayesplot)
# draws <- read_rds("output/k13snp_C469Y/bb_gne/draws.rds")
# mcmc_hist(draws)
library(cowplot)

kelch_pal <- c("#14B1E7", viridis(5))
partner_pal <- c("#c7047c", "#7ECE7E", "#174D97", "#E37210")

# r_hats <- coda::gelman.diag(draws,
#                             autoburnin = FALSE,
#                             multivariate = FALSE)
# message(summary(r_hats$psrf))

parameter_lookup <- list("gneiting_len" = "Lengthscale",
                         "gneiting_tim" = "Timescale",
                         "gneiting_sd" = "Kernel SD", # rename
                         "white_sd" = "White noise SD", # rename
                         "beta[1,1]" = "Intercept",
                         "beta[2,1]" = "Beta: time", # check which is which
                         "beta[3,1]" = "Beta: PfPR",
                         "rho" = "Rho")

# how did I work this out ?????
chains_to_exclude <- list("k13_marcse" = 3) # this particular chain was a bit dodgy


alldraws <- lapply(names(nice_name_lookup_all), function(marker){
  draws <- read_rds(paste0("output/", marker, "/bb_gne/draws.rds"))
  
  if (marker %in% names(chains_to_exclude)){
    draws <- draws[-c(chains_to_exclude[[marker]])]
  }
  
  draws %>% 
    do.call(what = rbind) %>% 
    as.data.frame() %>%
    rename_with(~ unlist(parameter_lookup), all_of(names(parameter_lookup))) %>%
    mutate(marker = nice_name_lookup_all[[marker]])
}) %>%
  do.call(what = bind_rows) %>%
  pivot_longer(cols = -c(marker), names_to = "var", values_to = "val") %>%
  mutate(marker = factor(marker, levels = nice_name_lookup_all))

set.seed(13517)

priors <- list(gneiting_len = sort(rnorm(4000, sd = 3)),
               gneiting_tim = sort(rnorm(4000, sd = 3)),
               gneiting_sd = sort(rnorm(4000, sd = 3)),
               white_sd = sort(rnorm(4000, sd = 3)),
               rho = rlnorm(4000),
               `beta[1,1]` = rnorm(4000),
               `beta[2,1]` = rnorm(4000),  
               `beta[3,1]` = rnorm(4000)) %>%
  do.call(what = bind_cols) %>%
  # sort out truncated normals jankily
  mutate(across(gneiting_len:white_sd, ~ ifelse(.x >= 0, .x, NA))) %>%
  drop_na() %>%
  rename_with(~ unlist(parameter_lookup), all_of(names(parameter_lookup))) %>%
  pivot_longer(everything(), values_to = "val", names_to = "var") %>%
  mutate(prio = "Priors")

nrow(priors)

# hyperparameters from model code:
# gneiting_len <- normal(0, 3, truncation = c(0, Inf))
# gneiting_tim <- normal(0, 3, truncation = c(0, Inf))
# gneiting_sd <- normal(0, 2, truncation = c(0, Inf))
# white_sd <- normal(0, 3, truncation = c(0, Inf)) # Median :1.041  Mean   :1.115
# rho <- greta::lognormal(0, 1)
# beta <- normal(0, 1, dim = 3)

library(ggh4x)

# Define individual limits for specific facet values

by_panel_x <- list(
  scale_x_continuous(limits = c(-1.5, 1.5)),
  scale_x_continuous(limits = c(-2, 4)),
  scale_x_continuous(limits = c(-4, 4)),
  scale_x_continuous(limits = c(0, 10)),
  scale_x_continuous(limits = c(0, 1.75)),
  scale_x_continuous(limits = c(0, 0.4)),
  scale_x_continuous(limits = c(0, 7)),
  scale_x_continuous(limits = c(0, 3.5))
)

p1 <- ggplot(data = alldraws %>% 
         filter(str_detect(marker, "Kelch"))) +
  geom_density(aes(x = val, col = marker), key_glyph = draw_key_path) +
  geom_density(data = priors, aes(x = val, linetype = prio), 
               key_glyph = draw_key_path) +
  facet_wrap(~ var, scales = "free") +
  scale_colour_manual(values = kelch_pal, "Kelch 13 markers") +
  scale_linetype_manual("", values = "dashed", guide = "none") +
  xlab("Parameter value") +
  ylab("Density") +
  facetted_pos_scales(
    x = by_panel_x
  ) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.85, 0.07),
        legend.spacing.y = unit(-0.75, "cm"),
        legend.background = element_rect(fill = NA),
        plot.margin = margin(1, 1, 25, 1),
        axis.title=element_text(size=10),
        legend.title=element_text(size=10))

p1

by_panel_x <- list(
  scale_x_continuous(limits = c(-1, 1)),
  scale_x_continuous(limits = c(-2, 2)),
  scale_x_continuous(limits = c(-4, 5)),
  scale_x_continuous(limits = c(0, 10)),
  scale_x_continuous(limits = c(0, 1.5)),
  scale_x_continuous(limits = c(0, 0.4)),
  scale_x_continuous(limits = c(0, 7)),
  scale_x_continuous(limits = c(0, 3.5))
)


# i would like to overlay posterior samples from multiple draws
# objects and prior samples
p2 <- ggplot(data = alldraws %>% 
               filter(marker %in% nice_name_lookup_main)) + # doesn't include aggregate
  geom_density(aes(x = val, col = marker), key_glyph = draw_key_path) +
  geom_density(data = priors, aes(x = val, linetype = prio), 
               key_glyph = draw_key_path) +
  facet_wrap(~ var, scales = "free") +
  scale_color_manual(values = partner_pal, "Partner drug markers") +
  xlab("Parameter value") +
  ylab("Density") +
  scale_linetype_manual("", values = "dashed") +
  facetted_pos_scales(
    x = by_panel_x
  ) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.85, 0.07),
        legend.spacing.y = unit(-0.5, "cm"),
        legend.background = element_rect(fill = NA),
        plot.margin = margin(1, 1, 25, 1),
        axis.title=element_text(size=10),
        legend.title=element_text(size=10)) +
  guides(colour = guide_legend(order = 1), 
           linetype = guide_legend(order = 2))

p2

p <- plot_grid(p1, p2, ncol = 1, align = "v")

ggsave("figures/posterior_densities.png", p, height = 10, width = 7.5)


