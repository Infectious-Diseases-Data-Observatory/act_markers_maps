# do some visualisation in here ...
setwd("~/Desktop/MARCSE/k13_seafrica")
source("code/setup.R")
source("code/vis/vis_funcs.R")
library(sf)

to_report <- "output/stats_to_report.txt"
reports <- "#### Vis script ####"

library(gridExtra)
library(grid)
library(cowplot)
library(patchwork)

# mask
afr <- world %>%
  filter(continent == "Africa") %>%
  vect() %>%
  crop(ext(-21, 63, -35, 37)) %>%
  st_as_sf()

#output_dir <- "output/circmat_pfmdr86/"

###############################################################################
# surveillance effort - very slow running atm ...

# p <- survey_effort_panel(c("output/k13/surveillance_effort_k13.grd",
#                       "output/crt76/surveillance_effort_crt76.grd",
#                       "output/mdr86/surveillance_effort_mdr86.grd"),
#                     lyr_names = c("Pfkelch13", "Pfcrt K76T", "Pfmdr1 N86Y"))
# 
# ggsave("figures/surveillance_effort_all.png", p,
#         height = 5, width=10, scale = 0.8)

p <- survey_effort_panel(c("output/k13_marcse/surveillance_effort_k13.grd",
                      "output/crt76/surveillance_effort_crt76.grd",
                      "output/mdr86/surveillance_effort_mdr86.grd"),
                    lyr_names = c("Kelch 13", "Pfcrt K76T", "Pfmdr1 N86Y"))

ggsave("figures/surveillance_effort_all_marcse.png", p, height = 6, width = 12)


# survey_effort_panel(c("output/k13/surveillance_effort_k13.grd",
#                       "output/k13_marcse/surveillance_effort_k13_marcse.grd"),
#                     lyr_names = c("moldm", "moldm + unpublished data"))
# ggsave("~/Desktop/presentations/MARCSE/surveillance_effort_marcse.png",
#        height = 5, width=8, scale = 0.9)

# after all that, I might just facet_wrap ....

###############################################################################
# time !

p1 <- pred_time_plot("output/k13_marcse/bb_gne/",
                     title = "(a) Kelch 13", show_pts = TRUE,
                     cut_year = 2028)
p2 <- pred_time_plot("output/crt76/bb_gne/",
                     title = "(b) Pfcrt K76T", show_pts = TRUE,
                     cut_year = 2028)
p3 <- pred_time_plot("output/mdr86/bb_gne/",
                     title = "(c) Pfmdr1 N86Y", show_pts = TRUE,
                     cut_year = 2028)
p4 <- pred_time_plot("output/mdr184/bb_gne/",
                     title = "(d) Pfmdr1 Y184F", show_pts = TRUE,
                     cut_year = 2028)
p5 <- pred_time_plot("output/mdr1246/bb_gne/",
                     title = "(e) Pfmdr1 D1246Y", show_pts = TRUE,
                     cut_year = 2028)

p1 + p2 + p3 + p4 + p5 +
  plot_layout(ncol = 1, guides = "collect", axis_title = "collect")
ggsave("figures/all_markers_time_bb.png", scale = 1.5, height = 7, width = 6)


p1 <- pred_time_plot_crint("output/k13_marcse/bb_gne/",
                     title = "(a) Kelch 13", show_pts = TRUE,
                     cut_year = 2028)
p2 <- pred_time_plot_crint("output/crt76/bb_gne/",
                     title = "(b) Pfcrt K76T", show_pts = TRUE,
                     cut_year = 2028)
p3 <- pred_time_plot_crint("output/mdr86/bb_gne/",
                     title = "(c) Pfmdr1 N86Y", show_pts = TRUE,
                     cut_year = 2028)
p4 <- pred_time_plot_crint("output/mdr184/bb_gne/",
                     title = "(d) Pfmdr1 Y184F", show_pts = TRUE,
                     cut_year = 2028)
p5 <- pred_time_plot_crint("output/mdr1246/bb_gne/",
                     title = "(e) Pfmdr1 D1246Y", show_pts = TRUE,
                     cut_year = 2028)
p1 + p2 + p3 + p4 + p5 +
  plot_layout(ncol = 1, guides = "collect", axis_title = "collect")
ggsave("figures/all_markers_time_bb_crint.png", scale = 1.5, height = 7, width = 6)

p1 <- pred_time_plot("output/k13_marcse/bb_gne/",
                     title = "(a) Pfkelch13",
                     cut_year = 2028)
p2 <- pred_time_plot("output/crt76/bb_gne/",
                     title = "(b) Pfcrt K76T",
                     cut_year = 2028)
p3 <- pred_time_plot("output/mdr86/bb_gne/",
                     title = "(c) Pfmdr1 N86Y",
                     cut_year = 2028)
p4 <- pred_time_plot("output/mdr184/bb_gne/",
                     title = "(d) Pfmdr1 Y184F",
                     cut_year = 2028)
p5 <- pred_time_plot("output/mdr1246/bb_gne/",
                     title = "(e) Pfmdr1 D1246Y",
                     cut_year = 2028)

p1 + p2 + p3 + p4 + p5 +
  plot_layout(ncol = 1, guides = "collect", axis_title = "collect")
ggsave("figures/all_markers_time_bb_nopts_astmh.png", scale = 1.5, height = 7, width = 6)

# make a version of this weighted by malaria case incidence
incid <- rast("data/incid_rasters_afr.tif")
names(incid) <- paste0("incid_", str_extract(names(incid), ".{4}$"))

# p1 <- pred_time_plot("output/k13_marcse/bb_gne/",
#                      title = "(a) Pfkelch13", incid = incid, ylab = "Annual cases")
# p2 <- pred_time_plot("output/crt76/bb_gne/",
#                      title = "(b) Pfcrt K76T", incid = incid, ylab = "Annual cases")
# p3 <- pred_time_plot("output/mdr86/bb_gne/",
#                      title = "(c) Pfmdr1 N86Y", incid = incid, ylab = "Annual cases")
# p4 <- pred_time_plot("output/mdr184/bb_gne/",
#                      title = "(d) Pfmdr1 Y184F", incid = incid, ylab = "Annual cases")
# p5 <- pred_time_plot("output/mdr1246/bb_gne/",
#                      title = "(e) Pfmdr1 D1246Y", incid = incid, ylab = "Annual cases")
# 
# p1 + p2 + p3 + p4 + p5 +
#   plot_layout(ncol = 1, guides = "collect", axis_title = "collect")
# ggsave("figures/all_markers_time_bb_incidence_weighted.png", scale = 1.5, height = 7, width = 6)
# ultimately this isn't super useful: e.g., for Kelch, it's saying that increasingly,
# the mutation intercepts with high transmission areas, but there are lots of pixels
# where incidence is very low ... really need to show quantiles of prevalence * annual incidence
# surface which will require a fairly different function



# p1 <- pred_time_factoring_incidence("output/k13_marcse/bb_gne/",
#                      title = "(a) Pfkelch13", incid = incid, ylab = "Annual cases",
#                      cut_year = 2028)
# p2 <- pred_time_factoring_incidence("output/crt76/bb_gne/",
#                      title = "(b) Pfcrt K76T", incid = incid, ylab = "Annual cases",
#                      cut_year = 2028)
# p3 <- pred_time_factoring_incidence("output/mdr86/bb_gne/",
#                      title = "(c) Pfmdr1 N86Y", incid = incid, ylab = "Annual cases",
#                      cut_year = 2028)
# p4 <- pred_time_factoring_incidence("output/mdr184/bb_gne/",
#                      title = "(d) Pfmdr1 Y184F", incid = incid, ylab = "Annual cases",
#                      cut_year = 2028)
# p5 <- pred_time_factoring_incidence("output/mdr1246/bb_gne/",
#                      title = "(e) Pfmdr1 D1246Y", incid = incid, ylab = "Annual cases",
#                      cut_year = 2028)
# 
# p1 + p2 + p3 + p4 + p5 +
#   plot_layout(ncol = 1, guides = "collect", axis_title = "collect")
# ggsave("figures/all_markers_time_bb_incidence_weighted.png", scale = 1.5, height = 7, width = 6)


# the aggregation factor here bothers me: estimate would probably be lower if I aggregated by median
# but would be losing tail ...
# should have slurmed this
# 2026 summary stats: 0.0584222927574479; (0.0294689492866093, 0.109705518437546)
# which is pretty close to what I reported with aggregation:

p1 <- pred_time_factoring_incidence("output/k13_marcse/bb_gne/",
                               title = "(a) Pfkelch13", incid = incid,
                               ylab = "Proportion of annual cases",
                               proportional = TRUE,
                               cut_year = 2028)
p2 <- pred_time_factoring_incidence("output/crt76/bb_gne/",
                               title = "(b) Pfcrt K76T", incid = incid,
                               ylab = "Proportion of annual cases",
                               proportional = TRUE, ylim = c(0,1),
                               cut_year = 2028)
p3 <- pred_time_factoring_incidence("output/mdr86/bb_gne/",
                               title = "(c) Pfmdr1 N86Y", incid = incid,
                               ylab = "Proportion of annual cases",
                               proportional = TRUE, ylim = c(0,1),
                               cut_year = 2028)
p4 <- pred_time_factoring_incidence("output/mdr184/bb_gne/",
                               title = "(d) Pfmdr1 Y184F", incid = incid,
                               ylab = "Proportion of annual cases",
                               proportional = TRUE, ylim = c(0,1),
                               cut_year = 2028)
p5 <- pred_time_factoring_incidence("output/mdr1246/bb_gne/",
                               title = "(e) Pfmdr1 D1246Y", incid = incid,
                               ylab = "Proportion of annual cases",
                               proportional = TRUE, ylim = c(0,1),
                               cut_year = 2028)

p1 + p2 + p3 + p4 + p5 +
  plot_layout(ncol = 1, guides = "collect", axis_title = "collect")
ggsave("figures/all_markers_time_bb_incidence_weighted_proportional.png", scale = 1.5, height = 7, width = 6)


plot(incid$incid_2024)
plot(sqrt(incid$incid_2024))
plot(log10(incid$incid_2024))
preds <- rast("output/k13_marcse/bb_gne/preds_medians.tif")
plot(preds$`2026_50`)
plot(sqrt(preds$`2026_50` * incid$incid_2024))
lower <- rast("output/k13_marcse/bb_gne/preds_lower.tif")
upper <- rast("output/k13_marcse/bb_gne/preds_upper.tif")

reports <- c(reports,
             paste("K13-affected cases -",
                   pred_time_factoring_incidence("output/k13_marcse/bb_gne/",
                                                 incid = incid,
                                                 proportional = TRUE,
                                                 cut_year = 2028, agg_fact = 1,
                                                 plot = FALSE)))

# make a version of this with highlights to k13 panel, highlights to other pred panels?
# after all that, not sure this adds much?
# zambezi <- list(xmin = 19, xmax = 26, ymin = -21, ymax = -14)
# victoria <- list(xmin = 27, xmax = 35, ymin = -6, ymax = 2)
# eswatini <- list(xmin = 34, xmax = 41, ymin = 12, ymax = 19)
#
# zoom_df <- rbind(zambezi, victoria, eswatini) %>%
#   as.data.frame() %>%
#   unnest() %>%
#   suppressWarnings()

# ended up abandoning this
# p1 <- pred_time_plot("output/k13/gneiting_sparse/preds_all.grd",
#                      title = "(a) Pfkelch13",
#                      zooms = zoom_df, zoom_pal = case_pal,
#                      points_path = "data/clean/moldm_k13_nomarker.csv")
# p2 <- pred_time_plot("output/crt76/gneiting_sparse/preds_all.grd",
#                      title = "(b) Pfcrt K76T",
#                      zooms = zoom_df, zoom_pal = case_pal,
#                      points_path = "data/clean/moldm_crt76.csv")
# p3 <- pred_time_plot("output/mdr86/gneiting_sparse/preds_all.grd",
#                      title = "(c) Pfmdr1 N86Y",
#                      zooms = zoom_df, zoom_pal = case_pal,
#                      points_path = "../moldm/clean/pfmdr_single_86.csv")
# p4 <- pred_time_plot("output/mdr184/gneiting_sparse/preds_all.grd",
#                      title = "(d) Pfmdr1 Y184F",
#                      zooms = zoom_df, zoom_pal = case_pal,
#                      points_path = "../moldm/clean/pfmdr_single_184.csv")
# p5 <- pred_time_plot("output/mdr1246/gneiting_sparse/preds_all.grd",
#                      title = "(e) Pfmdr1 D1246Y",
#                      zooms = zoom_df, zoom_pal = case_pal,
#                      points_path = "../moldm/clean/pfmdr_single_1246.csv")
#
# p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 1, guides = "collect", axis_title = "collect")
# ggsave("figures/all_markers_time_pts_zoom_gneiting.png", scale = 1.5, height = 7, width = 6)

#########################################################################
# wrap up model predictions for pfcrt/pfmdr1:

# executive decision: too many maps in the paper:
# putting uncertainties in supp .. also the code I wrote originally is too cumbersome

markers = c("crt76", "mdr86", "mdr184", "mdr1246")

years_to_plot <- as.character(seq(2002, 2026, 6))

to_plot <- lapply(markers, function(x){
  preds <- rast(paste0("output/", x, "/bb_gne/preds_medians.tif"))
  preds <- preds[[str_extract(names(preds), "^.{4}") %in% years_to_plot]]
  names(preds) <- paste0(str_extract(names(preds), "^.{4}"), "_", x)
  preds
}) %>%
  rast()

coords <- xyFromCell(to_plot, cells(to_plot))
vals <- terra::extract(to_plot, coords)
df <- cbind(coords, vals) %>%
  pivot_longer(starts_with("2"),
               names_to = "lyr",
               values_to = "val") %>%
  mutate(year = substr(lyr, 1, 4),
         marker = str_extract(lyr, "(?<=_).*"))

reports <- c(reports,
             paste("crt median prevalence 2026:",
                   median(values(to_plot$`2026_crt76`), na.rm=TRUE)))

p <- ggplot(df %>%
         mutate(marker = case_when(marker == "mdr1246"~ "Pfmdr1 D1246Y",
                                   marker == "mdr184" ~"Pfmdr1 Y184F",
                                   marker == "mdr86" ~ "Pfmdr1 N86Y",
                                   marker == "crt76" ~ "Pfcrt K76T"),
                marker = factor(marker,
                                levels = rev(c("Pfmdr1 D1246Y", "Pfmdr1 Y184F",
                                               "Pfmdr1 N86Y", "Pfcrt K76T"))))) +
  geom_tile(aes(x = x, y = y, fill = val)) +
  geom_sf(data = afr, fill = NA, linewidth = 0.1) +
  facet_grid(year ~ marker) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_gradientn(name = "Prevalence",
                       colors = blrd,
                       breaks = c(0, 0.5, 1),
                       labels = c("0\n(all wildtype)", "0.5", "1\n(all mutant)"),
                       limits = c(0,1)) +
  scale_x_continuous(breaks = seq(0, 40, 20), "Longitude") +
  scale_y_continuous(breaks = seq(-20, 40, 20), "Latitude") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(hjust = 0.5),
        legend.spacing.x = unit(4, "lines"),
        panel.spacing = unit(0, "lines")) +
  guides(fill = guide_colourbar(title.position = "top"),
         size = guide_legend(title.position = "top"))
ggsave("figures/crt_mdr_out_bb.png", p, height = 8, width = 5.5, scale = 1.1)

turbo_manip <- turbo(10)
turbo_manip <- c(turbo_manip[1:2],
  colorRampPalette(turbo_manip[3:8])(30),
  turbo_manip[9:10])



p <- ggplot(df %>%
              mutate(marker = case_when(marker == "mdr1246"~ "Pfmdr1 D1246Y",
                                        marker == "mdr184" ~"Pfmdr1 Y184F",
                                        marker == "mdr86" ~ "Pfmdr1 N86Y",
                                        marker == "crt76" ~ "Pfcrt K76T"),
                     marker = factor(marker,
                                     levels = rev(c("Pfmdr1 D1246Y", "Pfmdr1 Y184F",
                                                    "Pfmdr1 N86Y", "Pfcrt K76T"))))) +
  geom_tile(aes(x = x, y = y, fill = val)) +
  geom_sf(data = afr, fill = NA, linewidth = 0.1) +
  facet_grid(year ~ marker) +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_gradientn(name = "Prevalence",
                       breaks = c(0, 0.5, 1),
                       labels = c("0\n(all wildtype)", "0.5", "1\n(all mutant)"),
                       limits = c(0,1),
                       colours = turbo_manip) +
  scale_x_continuous(breaks = seq(0, 40, 20), "Longitude") +
  scale_y_continuous(breaks = seq(-20, 40, 20), "Latitude") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(hjust = 0.5),
        legend.spacing.x = unit(4, "lines"),
        panel.spacing = unit(0, "lines")) +
  guides(fill = guide_colourbar(title.position = "top"),
         size = guide_legend(title.position = "top"))
ggsave("figures/crt_mdr_out_bb_turbo.png", p, height = 8, width = 5.5, scale = 1.1)

# horizontal for presentations folder
p <- ggplot(df %>%
              # filter(year %in% c("2002", "2014", "2026")) %>%
              mutate(marker = case_when(marker == "mdr1246"~ "Pfmdr1 D1246Y",
                                        marker == "mdr184" ~"Pfmdr1 Y184F",
                                        marker == "mdr86" ~ "Pfmdr1 N86Y",
                                        marker == "crt76" ~ "Pfcrt K76T"),
                     marker = factor(marker,
                                     levels = rev(c("Pfmdr1 D1246Y", "Pfmdr1 Y184F",
                                                    "Pfmdr1 N86Y", "Pfcrt K76T"))))) +
  geom_tile(aes(x = x, y = y, fill = val)) +
  geom_sf(data = afr, fill = NA, linewidth = 0.1) +
  facet_grid(marker~year, switch = "y") +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_gradientn(name = "Prevalence",
                       colors = blrd,
                       breaks = c(0, 0.5, 1),
                       labels = c("0 (all wildtype)", "0.5", "1 (all mutant)"),
                       limits = c(0,1)) +
  scale_x_continuous(breaks = seq(0, 40, 20), "Longitude") +
  scale_y_continuous(breaks = seq(-20, 40, 20), "Latitude") +
  theme_bw() +
  theme(#legend.position = "bottom",
        #legend.title = element_text(hjust = 0.5),
        #legend.spacing.x = unit(4, "lines"),
        panel.spacing = unit(0, "lines"))
ggsave("~/Desktop/presentations/marcse/crt_mdr_out_bb.png", p, 
       height = 5, width = 7, scale = 1.2)

#########################################################################
# here is uncertainties:
to_plot <- lapply(markers, function(x){
  preds <- rast(paste0("output/", x, "/bb_gne/preds_sds.tif"))
  preds <- preds[[str_extract(names(preds), "^.{4}") %in% years_to_plot]]
  names(preds) <- paste0(str_extract(names(preds), "^.{4}"), "_", x)
  preds
}) %>%
  rast()

coords <- xyFromCell(to_plot, cells(to_plot))
vals <- terra::extract(to_plot, coords)
df <- cbind(coords, vals) %>%
  pivot_longer(starts_with("2"),
               names_to = "lyr",
               values_to = "val") %>%
  mutate(year = substr(lyr, 1, 4),
         marker = str_extract(lyr, "(?<=_).*")) %>%
  mutate(marker = case_when(marker == "mdr1246"~ "Pfmdr1 D1246Y",
                            marker == "mdr184" ~"Pfmdr1 Y184F",
                            marker == "mdr86" ~ "Pfmdr1 N86Y",
                            marker == "crt76" ~ "Pfcrt K76T"),
         marker = factor(marker,
                         levels = rev(c("Pfmdr1 D1246Y", "Pfmdr1 Y184F",
                                        "Pfmdr1 N86Y", "Pfcrt K76T"))))


p1 <- sd_panelled_plot(df, "Pfcrt K76T", lab = "Estimate SD") +
  theme(strip.text.y = element_blank(),
        strip.background.y = element_blank(),
        axis.title.x = element_blank())
p2 <- sd_panelled_plot(df, "Pfmdr1 N86Y", lab = "Estimate SD") +
  theme(strip.text.y = element_blank(),
        strip.background.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank())
p3 <- sd_panelled_plot(df, "Pfmdr1 Y184F", lab = "Estimate SD") +
  theme(strip.text.y = element_blank(),
        strip.background.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank())
p4 <- sd_panelled_plot(df, "Pfmdr1 D1246Y", lab = "Estimate SD") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank())

# p <- ggplot(df %>%
#               mutate(marker = case_when(marker == "mdr1246"~ "Pfmdr1 D1246Y",
#                                         marker == "mdr184" ~"Pfmdr1 Y184F",
#                                         marker == "mdr86" ~ "Pfmdr1 N86Y",
#                                         marker == "crt76" ~ "Pfcrt K76T"),
#                      marker = factor(marker,
#                                      levels = rev(c("Pfmdr1 D1246Y", "Pfmdr1 Y184F",
#                                                     "Pfmdr1 N86Y", "Pfcrt K76T"))))) +
#   geom_tile(aes(x = x, y = y, fill = val)) +
#   geom_sf(data = afr, fill = NA, linewidth = 0.1) +
#   facet_grid(year ~ marker) +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   scale_fill_gradientn(name = "Estimate SD",
#                        colors = oranges,
#                        breaks = seq(0.05, 0.25, length.out = 3)) + # they got a bit squished so manually setting
#   scale_x_continuous(breaks = seq(-20, 40, 20), "Longitude") +
#   scale_y_continuous(breaks = seq(-20, 40, 20), "Latitude") +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         legend.title = element_text(hjust = 0.5),
#         legend.spacing.x = unit(4, "lines"),
#         panel.spacing = unit(0, "lines")) +
#   guides(fill = guide_colourbar(title.position = "top"),
#          size = guide_legend(title.position = "top"))
# ggsave("figures/crt_mdr_out_bb_sds.png", p, height = 9, width = 6.5)


text_patch <- ggplot() + 
  annotate("text", x = 1, y = 1, label = "Longitude", size = 4, fontface = "plain") +
  theme_void() +
  coord_cartesian(clip = "off")

p <- (p1 | p2 | p3 | p4) +
  inset_element(text_patch, left = -2.5, bottom = -0.05, right = 0, top = 0.05)


ggsave("figures/crt_mdr_out_bb_sds.png", p, height = 9.5, width = 7)


#########################################################################
# and unscaled uncertainties
to_plot <- lapply(markers, function(x){
  preds <- rast(paste0("output/", x, "/bb_gne/preds_sdscaled.tif"))
  preds <- preds[[str_extract(names(preds), "^.{4}") %in% years_to_plot]]
  names(preds) <- paste0(str_extract(names(preds), "^.{4}"), "_", x)
  preds
}) %>%
  rast()

coords <- xyFromCell(to_plot, cells(to_plot))
vals <- terra::extract(to_plot, coords)
df <- cbind(coords, vals) %>%
  pivot_longer(starts_with("2"),
               names_to = "lyr",
               values_to = "val") %>%
  mutate(year = substr(lyr, 1, 4),
         marker = str_extract(lyr, "(?<=_).*")) %>%
  mutate(marker = case_when(marker == "mdr1246"~ "Pfmdr1 D1246Y",
                            marker == "mdr184" ~"Pfmdr1 Y184F",
                            marker == "mdr86" ~ "Pfmdr1 N86Y",
                            marker == "crt76" ~ "Pfcrt K76T"),
         marker = factor(marker,
                         levels = rev(c("Pfmdr1 D1246Y", "Pfmdr1 Y184F",
                                        "Pfmdr1 N86Y", "Pfcrt K76T"))))

sd_panelled_plot <- function(df, mark, lab = "Estimate SD\n(unscaled)"){
  # giving each model its own colour bar: split out into separate plots
  # and plot_grid back together
  ggplot(df %>%
           filter(marker == mark)) +
    geom_tile(aes(x = x, y = y, fill = val)) +
    geom_sf(data = afr, fill = NA, linewidth = 0.1) +
    facet_grid(year ~ marker) +
    xlab("Longitude") +
    ylab("Latitude") +
    scale_fill_gradientn(name = lab,
                         colors = oranges) + # they got a bit squished so manually setting
    scale_x_continuous(breaks = seq(-20, 40, 20), "Longitude") +
    scale_y_continuous(breaks = seq(-20, 40, 20), "Latitude") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_text(hjust = 0.5),
          legend.spacing.x = unit(4, "lines"),
          panel.spacing = unit(0, "lines")) +
    guides(fill = guide_colourbar(title.position = "top"),
           size = guide_legend(title.position = "top"))
}

p1 <- sd_panelled_plot(df, "Pfcrt K76T") +
  theme(strip.text.y = element_blank(),
        strip.background.y = element_blank(),
        axis.title.x = element_blank())
p2 <- sd_panelled_plot(df, "Pfmdr1 N86Y") +
  theme(strip.text.y = element_blank(),
        strip.background.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank())
p3 <- sd_panelled_plot(df, "Pfmdr1 Y184F") +
  theme(strip.text.y = element_blank(),
        strip.background.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank())
p4 <- sd_panelled_plot(df, "Pfmdr1 D1246Y") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank())

# p <- ggplot(df) +
#   geom_tile(aes(x = x, y = y, fill = val)) +
#   geom_sf(data = afr, fill = NA, linewidth = 0.1) +
#   facet_grid(year ~ marker) +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   scale_fill_gradientn(name = "Estimate SD (unscaled)",
#                        colors = oranges) + # they got a bit squished so manually setting
#   scale_x_continuous(breaks = seq(-20, 40, 20), "Longitude") +
#   scale_y_continuous(breaks = seq(-20, 40, 20), "Latitude") +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         legend.title = element_text(hjust = 0.5),
#         legend.spacing.x = unit(4, "lines"),
#         panel.spacing = unit(0, "lines")) +
#   guides(fill = guide_colourbar(title.position = "top"),
#          size = guide_legend(title.position = "top"))

# p <- plot_grid(p1, p2, p3, p4, nrow = 1, align = "v")
               #rel_widths = c(1.4, 1, 1, 1.2))

library(patchwork)

text_patch <- ggplot() + 
  annotate("text", x = 1, y = 1, label = "Longitude", size = 4, fontface = "plain") +
  theme_void() +
  coord_cartesian(clip = "off")

p <- (p1 | p2 | p3 | p4) +
  inset_element(text_patch, left = -2.5, bottom = -0.05, right = 0, top = 0.05)


ggsave("figures/crt_mdr_out_bb_sdscaled.png", p, height = 9.5, width = 7)

# p1 <- map_pred_row("output/crt76/bb_gne/preds_medians.tif",
#                    exte = afr,
#                    years = years_to_plot, field = "50", pal = blrd,
#                    legend_lim = c(0,1), ylab = "", top_pan = TRUE)#
# p2 <- map_pred_row("output/mdr86/bb_gne/preds_medians.tif",
#                    years = years_to_plot, field = "50", pal = blrd,
#                    legend_lim = c(0,1), xlab = "", ylab = "")# ylab = "(a)")
# p3 <- map_pred_row("output/mdr184/bb_gne/preds_medians.tif",
#                    years = years_to_plot, field = "50", pal = blrd,
#                    legend_lim = c(0,1), xlab = "", ylab = "")#  ylab = "(c)")
# p4 <- map_pred_row("output/mdr1246/bb_gne/preds_medians.tif",
#                    years = years_to_plot, field = "50", pal = blrd,
#                    legend_lim = c(0,1), xlab = "",ylab = "")#  ylab = "(e)")
# p5 <- map_pred_row("output/mdr86/bb_gne/preds_sds.tif",
#                    years = years_to_plot, field = "sd", pal = oranges,
#                    legend_lim = c(0,0.2), xlab = "",ylab = "")#  ylab = "(b)")
# p6 <- map_pred_row("output/mdr86/bb_gne/preds_sds.tif",
#                    years = years_to_plot, field = "sd", pal = oranges,
#                    legend_lim = c(0,0.2), xlab = "",ylab = "")#  ylab = "(b)")
# p7 <- map_pred_row("output/mdr184/bb_gne/preds_sds.tif",
#                    years = years_to_plot, field = "sd", pal = oranges,
#                    legend_lim = c(0,0.2), xlab = "",ylab = "")#  ylab = "(d)")
# p8 <- map_pred_row("output/mdr1246/bb_gne/preds_sds.tif",
#                    years = years_to_plot, field = "sd", pal = oranges,
#                    legend_lim = c(0,0.2), xlab = "",ylab = "")#  ylab = "(f)")
# # legend_lim was c(0, 0.2) for the main fig ... although this probably could have been less



# # chuck the rows together
# pcol <- plot_grid(p2 + theme(legend.position = "none"),
#              p6 + theme(legend.position = "none"),
#              p3 + theme(legend.position = "none"),
#              p7 + theme(legend.position = "none"),
#              p4 + theme(legend.position = "none"),
#              p8 + theme(legend.position = "none"),
#              #ncol = 1, rel_heights = c(1.21, rep(1, 5))) +
#              ncol = 1, rel_heights = c(1.167, rep(1, 5))) +
#   theme(panel.spacing = unit(0, "cm"))
#
# #pcol
#
# #pcol <- grid.arrange(arrangeGrob(pcol, left = y.grob, bottom = x.grob))
#
# legend1 = get_legend(p2)
# legend2 = get_legend(p6)
#
# plegend <- plot_grid(legend1, legend2, ncol = 1, axis = "r")
#
# #plegend
#
# p <- plot_grid(pcol, plegend, rel_widths = c(0.8,0.2))
#
# # this df worked when I was using subfigure labels ...
# # df <- data.frame(xmin = rep(0.02, 3),
# #                  xmax = rep(0.062, 3),
# #                  ymin = c(0.017, 0.34, 0.66),
# #                  ymax = c(0.312, 0.635, 0.955),
# #                  lab = c("D1246Y", "Y184F", "N86Y"))
# df <- data.frame(xmin = rep(0.001, 3),
#                  xmax = rep(0.028, 3),
#                  ymin = c(0.0155, 0.337, 0.663),
#                  ymax = c(0.315, 0.639, 0.965),
#                  lab = c("D1246Y", "Y184F", "N86Y"))
# p <- p +
#   # tried adding outer margin but that didn't do anything
#   geom_rect(data = df, aes(xmin=xmin, xmax=xmax, ymin=ymin,
#                           ymax=ymax),
#             colour="grey10", fill="grey85", linewidth=0.3) +
#   geom_text(data = df, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
#                             label = lab), angle = 90) +
#   geom_text(data = data.frame(x = rep(0.85, 2), y = c(0.34, 0.84),
#                               label = c("Estimate SD", "Prevalence")),
#             aes(x = x, y = y, label = label))
#
# p
#
# # gave up on add_sub
# ggsave("figures/mdr_out_bb.png", height = 9, width = 8.5)

########################################################################
# can I fit crt in too?
# p2 <- map_pred_row("output/mdr86/bb_gne/preds_all.tif",
#                    years = years_to_plot, field = "50", pal = blrd,
#                    legend_lim = c(0,1), xlab = "", ylab = "")# ylab = "(a)")

# chuck the rows together
# pcol <- plot_grid(p1 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
#                   p5 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
#                   p2 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
#                   p6 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
#                   p3 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
#                   p7 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
#                   p4 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
#                   p8 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
#                   #ncol = 1, rel_heights = c(1.21, rep(1, 5))) +
#                   ncol = 1, rel_heights = c(1.22, rep(1, 7))) +
#   theme(panel.spacing = unit(0, "cm"))
#
# #pcol
#
# #pcol <- grid.arrange(arrangeGrob(pcol, left = y.grob, bottom = x.grob))
#
# legend1 = get_legend(p2)
# legend2 = get_legend(p6)
#
# plegend <- plot_grid(legend1, legend2, ncol = 1, axis = "r")
#
# #plegend
#
# p <- plot_grid(pcol, plegend, rel_widths = c(0.8,0.22))
#
# # this df worked when I was using subfigure labels ...
# # df <- data.frame(xmin = rep(0.02, 3),
# #                  xmax = rep(0.062, 3),
# #                  ymin = c(0.017, 0.34, 0.66),
# #                  ymax = c(0.312, 0.635, 0.955),
# #                  lab = c("D1246Y", "Y184F", "N86Y"))
# df <- data.frame(xmin = rep(0.001, 4),
#                  xmax = rep(0.035, 4),
#                  ymin = c(0.005, 0.249, 0.492, 0.735),
#                  ymax = c(0.242, 0.487, 0.73, 0.975),
#                  lab = c("Pfmdr1 D1246Y", "Pfmdr1 Y184F", "Pfmdr1 N86Y", "Pfcrt K76T"))
# p <- p +
#   # tried adding outer margin but that didn't do anything
#   geom_rect(data = df, aes(xmin=xmin, xmax=xmax, ymin=ymin,
#                            ymax=ymax),
#             colour="grey10", fill="grey85", linewidth=0.3) +
#   geom_text(data = df, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
#                            label = lab), angle = 90) +
#   geom_text(data = data.frame(x = rep(0.85, 2), y = c(0.34, 0.84),
#                               label = c("Estimate SD", "Prevalence")),
#             aes(x = x, y = y, label = label))
#
# ggsave("figures/crt_mdr_out_bb.png", p, height = 9, width = 7)

##############################################################################
# need a short version for powerpoint

# chuck the rows together
# pcol <- plot_grid(p1 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
#                   p2 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
#                   p3 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
#                   p4 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
#                   ncol = 1, rel_heights = c(1.22, rep(1, 3))) +
#   theme(panel.spacing = unit(0, "cm"))

# years_to_plot <- c(2006, 2014, 2022)
# preds <- rast(paste0("output/", c("crt76", "mdr86", "mdr184", "mdr1246"), "/gneiting_sparse/preds_all.grd")) %>%
#   aggregate(fact=3)
# preds <- preds[[grep("median", names(preds))]]
# names(preds) <- paste0(names(preds), "_", rep(c("crt76", "mdr86", "mdr184", "mdr1246"), each = 25))
# coords <- xyFromCell(preds, cells(preds))
# vals <- terra::extract(preds, coords)
# df <- cbind(coords, vals) %>%
#   pivot_longer(starts_with("2"),
#                names_to = "lyr",
#                values_to = "val") %>%
#   mutate(year = substr(lyr, 1, 4),
#          tag = substr(lyr, 11, 14),
#          marker = str_extract(lyr, "[^_]+$")) %>%
#   filter(year %in% years_to_plot & tag == "medi") # pick out year and thingo
#
# df <- df %>% mutate(marker = case_when(marker == "crt76" ~ "Pfcrt K76T",
#                                        marker == "mdr86" ~ "Pfmdr1 N86Y",
#                                        marker == "mdr184" ~ "Pfmdr1 Y184F",
#                                        marker == "mdr1246" ~ "Pfmdr1 D1246Y")) %>%
#   mutate(marker = factor(marker, levels = c("Pfcrt K76T", "Pfmdr1 N86Y", "Pfmdr1 Y184F", "Pfmdr1 D1246Y")))
#
# ggplot() +
#   geom_sf(data = st_as_sf(afr), fill = "white") +
#   geom_tile(aes(x = x, y = y, fill = val), data = df) +
#   scale_fill_gradientn(name = "Prevalence",
#                        colors = blrd,
#                        breaks = c(0, 0.5, 1),
#                        labels = c("0  (all wildtype)", "0.5", "1  (all mutant)"),
#                        limits = c(0,1)) +
#   facet_grid(year ~ marker) +
#   scale_x_continuous(breaks = seq(-20, 40, 20), "Longitude") +
#   scale_y_continuous(breaks = seq(-20, 40, 20), "Latitude") +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   theme_bw()
#
# ggsave("~/Desktop/presentations/MARCSE/crt_mdr_out.png", scale = 1.7, height = 3, width = 5)

# years_to_plot <- c(2006, 2014, 2022)
# preds <- rast(paste0("output/", c("crt76", "mdr86", "mdr184", "mdr1246"), "/gneiting_sparse/preds_all.grd")) %>%
#   aggregate(fact=3)
# preds <- preds[[grep("sd", names(preds))]]
# names(preds) <- paste0(names(preds), "_", rep(c("crt76", "mdr86", "mdr184", "mdr1246"), each = 25))
# coords <- xyFromCell(preds, cells(preds))
# vals <- terra::extract(preds, coords)
# df <- cbind(coords, vals) %>%
#   pivot_longer(starts_with("2"),
#                names_to = "lyr",
#                values_to = "val") %>%
#   mutate(year = substr(lyr, 1, 4),
#          tag = substr(lyr, 11, 14),
#          marker = str_extract(lyr, "[^_]+$"))
# df <- filter(df, year %in% years_to_plot) # pick out year and thingo
#
# df <- df %>% mutate(marker = case_when(marker == "crt76" ~ "Pfcrt K76T",
#                                        marker == "mdr86" ~ "Pfmdr1 N86Y",
#                                        marker == "mdr184" ~ "Pfmdr1 Y184F",
#                                        marker == "mdr1246" ~ "Pfmdr1 D1246Y")) %>%
#   mutate(marker = factor(marker, levels = c("Pfcrt K76T", "Pfmdr1 N86Y", "Pfmdr1 Y184F", "Pfmdr1 D1246Y")))
#
# ggplot() +
#   geom_sf(data = st_as_sf(afr), fill = "white") +
#   geom_tile(aes(x = x, y = y, fill = val), data = df) +
#   scale_fill_gradientn(name = "Uncertainty",
#                        colors = oranges, trans="sqrt") +
#   facet_grid(year ~ marker) +
#   scale_x_continuous(breaks = seq(-20, 40, 20), "Longitude") +
#   scale_y_continuous(breaks = seq(-20, 40, 20), "Latitude") +
#   theme_bw() +
#   xlab("Longitude") +
#   ylab("Latitude")
#
# ggsave("~/Desktop/presentations/MARCSE/crt_mdr_out_sd.png", scale = 1.7, height = 3, width = 5)

##############################################################################

# coords <- xyFromCell(pfpr, cells(pfpr))
# vals <- terra::extract(pfpr, coords)
# df <- cbind(coords, vals) %>%
#   mutate(cell = 1:nrow(coords)) %>%
#   pivot_longer(starts_with("pfpr"),
#                names_to = "lyr",
#                values_to = "val") %>%
#   mutate(year = as.numeric(substr(lyr, 6, 9))) # pick out year and thingo
#
# cells_to_plot <- sample(1:nrow(coords), 50, replace = FALSE)
#
# df_summ <- df %>%
#   group_by(year) %>%
#   summarise(q = list(quantile(val, c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)))) %>%
#   unnest_wider(q) %>%
#   ungroup() %>%
#   mutate(year = as.numeric(year))
#
# pal = iddoPal::iddo_palettes$soft_blues
#
# p <- ggplot(df_summ) +
#   geom_line(aes(x = year, y = `0%`, linetype = "0% - 100%")) +
#   geom_line(aes(x = year, y = `100%`, linetype = "0% - 100%")) +
#   geom_ribbon(aes(x = year, ymin = `2.5%`, ymax = `97.5%`, fill = "2.5% - 97.5%")) + #fill=pal[6]) +
#   geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`, fill = "25% - 75%")) + #fill=pal[4]) +
#   geom_ribbon(aes(x = year, ymin = `50%`, ymax = `50%`, fill = "50%")) + #fill=pal[1]) +
#   geom_line(aes(x = year, y = `50%`), col = pal[1], linewidth = 1) +
#   scale_linetype_manual("", values = c("0% - 100%" = 2)) +
#   scale_fill_manual("", values = c("2.5% - 97.5%" = pal[6], "25% - 75%" = pal[4], "50%" = pal[1])) +
#   geom_line(data = df %>% filter(cell %in% cells_to_plot),
#             mapping = aes(x = year, y = val, colour = cell, group = cell)) +
#   ylab("PfPR") +
#   xlab("Year") +
#   theme_bw() +
#   theme(legend.spacing.y = unit(-10, "cm"),
#         legend.background = element_rect(fill = NA))
#
# p

# P441L in SAfr
preds <- rast("output/k13snp_P441L/bb_gne/preds_medians.tif")
tmp <- mask(preds$`2026_50`, afr %>% filter(iso_a3 %in% c("ZMB", "ZWE", "NAM", "AGO", "BWA"))) %>%
  trim()
plot(tmp)

reports <- c(reports,
             "max P441L prevalence SAfr 2026: ", max(values(tmp), na.rm = TRUE))

# plot individual markers

dirs <- list.dirs("output")
snps <- dirs[grepl("k13snp", dirs) & str_count(dirs, "/") == 1]
snps <- str_extract(snps, "(?<=_).*")

# will need to open and trim to custom extents
years_to_plot <- as.character(seq(2014, 2026, 3))

exts <- list(A675V = c(22, 38, -3, 13), #c(20, 40, -15, 5),
             C469Y = c(22, 38, -3, 13),
             P441L = c(15, 38, -22, 5),
             R561H = c(22, 38, -13, 3),
             R622I = c(32, 42, 8, 18))

to_plot <- lapply(snps, function(snp){
  message(snp)
  preds <- rast(paste0("output/k13snp_", snp, "/bb_gne/preds_medians.tif"))
  preds <- preds[[str_extract(names(preds), "^.{4}") %in% years_to_plot]]
  preds <- crop(preds, exts[[snp]])

  names(preds) <- paste0(str_extract(names(preds), "^.{4}"), "_", snp)
  as.data.frame(preds, xy = TRUE) %>%
    mutate(snp = snp) %>%
    pivot_longer(starts_with("2"),
                 names_to = "lyr",
                 values_to = "val") %>%
    mutate(year = substr(lyr, 1, 4))
})

df = do.call(bind_rows, to_plot)

snp_subplot <- function(df,
                        marker,
                        strip.text = element_blank(),
                        strip.background = element_blank(),
                        axis.title.x = element_blank(),
                        xbreaks = waiver(),
                        ybreaks = waiver(),
                        fill_lims = c(0,0.657)){
  ggplot(data = df %>% filter(snp == marker)) +
    geom_tile(aes(x = x, y = y, fill = val)) +
    geom_sf(data = afr %>%
              st_crop(ext(exts[[marker]])), fill = NA, col = "grey80") +
    facet_grid(snp ~ year) +
    scale_fill_viridis_c("Prevalence", limits = fill_lims, trans = "sqrt",
                        na.value = "white") +
    scale_x_continuous(breaks = xbreaks, "Longitude") +
    scale_y_continuous(breaks = ybreaks) +
    theme_bw() +
    theme(strip.text.x = strip.text,
          strip.background.x = strip.background,
          axis.title.x = axis.title.x,
          axis.title.y = element_blank(),
          panel.spacing = unit(0, "lines"),
          plot.margin = margin(2, 2, 0, 0),
          legend.position = "bottom",
          legend.title = element_text(hjust = 0.5, vjust = 0.8)) +
    guides(fill = guide_colourbar(barwidth = unit(5, "cm")))
}

p1 <- snp_subplot(df, "A675V",
                  strip.text = element_text(colour = "grey10",
                                            size   = 11,
                                            margin = margin(5.5, 5.5, 5.5, 5.5)),
                  strip.background = element_rect(
                    fill   = "grey85",
                    colour = "grey20"
                  ),
                  xbreaks = seq(20, 35, 5)) %>%
  suppressWarnings()
p2 <- snp_subplot(df, "C469Y") %>%
  suppressWarnings()
p3 <- snp_subplot(df, "P441L",
                  xbreaks = seq(15, 35, 10)) %>%
  suppressWarnings()
p4 <- snp_subplot(df, "R561H") %>%
  suppressWarnings()
p5 <- snp_subplot(df, "R622I",
                  axis.title.x = element_text(),
                  xbreaks = seq(32, 40, 4)) %>%
  suppressWarnings()

p <- plot_grid(p1 + theme(legend.position = "none"),
          p2 + theme(legend.position = "none"),
          p3 + theme(legend.position = "none"),
          p4 + theme(legend.position = "none"),
          p5,
          # alas rel_widths = rep(1, 5), doesn't do the job
          rel_heights = c(1.21, 1, 1.2, 1, 1.58),
          ncol = 1)

pl <- ggdraw(p + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 15))) +
  draw_label("Latitude", x = 0.01, y = 0.5,
    hjust = 0.5, vjust = 0.5, angle = 90, size = 12
  )
ggsave("figures/k13_snps.png", pl, height = 7.55, width = 6, scale = 1.3)

# ggsave("figures/k13_snps.png", p, height = 7.75, width = 6, scale = 1.3)

ggplot(df %>% filter(year == "2026")) +
  geom_tile(aes(fill = val, x = x, y = y)) +
  scale_fill_viridis_c(trans = "sqrt", limits = c(0, 0.3)) +
  facet_wrap(~snp)

# just wondering if this might not be more informative as a contour plot ...

##############################################################################
# compare obs to pred
# could also plot obs against prediction uncertainty

#mut_data$pred[is.na(mut_data$pred)] = 0.4

cex_transform <- function(from){
  sqrt(from) / sqrt(max(from))
  #log10(from) / log10(max(from)) *4
}




##############################################################################
# the same but ggplot

# gosh I prefer how base plotting deals with panels :/
# manipulate these further when I decide whether I want facetting etc
# obs_prev_panel("data/clean/moldm_k13_nomarker.csv",
#                "output/k13/circmat_sparse/preds_all.grd",
#                "k13 circmat", xlim = c(0, 0.4), ylim = c(0, 0.4))#,
#                #facet_bins = c(2010, 2014, 2021))

# obs_prev_panel("data/clean/moldm_k13_nomarker.csv",
#                 "output/k13/gneiting_sparse/preds_all.grd",
#                 "k13 gneiting", xlim = c(0, 0.4), ylim = c(0, 0.4))

# obs_prev_panel("data/clean/moldm_marcse_k13_nomarker.csv",
#                "output/k13_marcse/gneiting_ahmc/preds_all.tif",
#                "", xlim = c(0, 0.6), ylim = c(0, 0.6),
#                ave_tag = "_50")
# ggsave("figures/residuals_k13m_gne.png", height = 3, width = 4, scale = 2)
# 
# obs_prev_panel("data/clean/pfmdr_single_86.csv",
#                 "output/mdr86/gneiting_ahmc/preds_all.tif",
#                 "mdr86 gneiting",
#                ave_tag = "_50") #, facet_bins = c(2008, 2012, 2016, 2020))
# ggsave("figures/residuals_86_gne.png", height = 9, width = 5, scale = 2)

# obs_prev_panel("data/clean/pfmdr_single_86.csv",
#                "output/mdr86/bb_gne/preds_all.tif",
#                "mdr86 gneiting",
#                ave_tag = "_50")
# ggsave("figures/residuals_86_bb.png", height = 9, width = 5, scale = 2)

# obs_prev_panel("data/clean/pfmdr_single_86.csv",
#                "output/mdr86/circmat/preds_all.grd", "mdr86 circmat",
#                ave_tag = "_post_median")
#                #facet_bins = c(2008, 2012, 2016, 2020) )
# ggsave("~/Desktop/residuals_86.png", height = 9, width = 5, scale = 2)

# obs_prev_panel("data/clean/pfmdr_single_1246.csv",
#                 "output/mdr1246/gneiting_ahmc/preds_all.tif",
#                #facet_bins = c(2008, 2012, 2016, 2020),
#                 "mdr1246 gneiting",
#                ave_tag = "_50")
# ggsave("figures/residuals_1246_gne.png", height = 9, width = 5, scale = 2)

# obs_prev_panel("data/clean/pfmdr_single_1246.csv",
#                "output/mdr1246/bb_gne/preds_all.tif",
#                #facet_bins = c(2008, 2012, 2016, 2020),
#                "mdr1246 gneiting", buffer = 100000)
# ggsave("figures/residuals_1246_bb.png", height = 9, width = 5, scale = 2)

# obs_prev_panel("data/clean/pfmdr_single_1246.csv",
#                "output/mdr1246/circmat/preds_all.grd",
#                #facet_bins = c(2008, 2012, 2016, 2020),
#                "mdr1246 circmat",
#                ave_tag = "_post_median")
# ggsave("~/Desktop/residuals_1246.png", height = 9, width = 5, scale = 2)

# obs_prev_panel("data/clean/pfmdr_single_184.csv",
#                 "output/mdr184/gneiting_ahmc/preds_all.tif",
#                 "mdr184 gneiting",
#                ave_tag = "_50")
# ggsave("figures/residuals_184_gne.png", height = 9, width = 5, scale = 2)

# oops got rid of this
# obs_prev_panel("data/clean/pfmdr_single_184.csv",
#                "output/mdr184/circmat/preds_all.grd",
#                 "mdr184 circmat",
#                ave_tag = "_post_median")

# obs_prev_panel("data/clean/pfmdr_single_184.csv",
#                "output/mdr184/bb_gne/preds_all.tif",
#                "mdr184 gneiting", buffer = 100000)
# ggsave("figures/residuals_184_bb.png", height = 9, width = 5, scale = 2)

# how do I capture clustered data where the model picks somewhere in the middle
# but there's lots of variance at the same location?
# rolling window?
# this is where the CV would be helpful to make my point ...

##############################################################################
# visualise coverages - relocating coverages.png to validation script

# here's bb_gne all models crints

p1 <- lower_upper_panel_crints("output/k13_marcse/bb_gne/")
p2 <- lower_upper_panel_crints("output/crt76/bb_gne/")
p3 <- lower_upper_panel_crints("output/mdr86/bb_gne/")
p4 <- lower_upper_panel_crints("output/mdr184/bb_gne/")
p5 <- lower_upper_panel_crints("output/mdr1246/bb_gne/")


p <- plot_grid(p1 + theme(legend.position = "none", axis.title = element_blank()),
               p2 + theme(legend.position = "none", axis.title = element_blank()),
               p3 + theme(legend.position = "none", axis.title = element_blank()),
               p4 + theme(legend.position = "none", axis.title = element_blank()),
               p5 + theme(legend.position = "none", axis.title = element_blank()),
               nrow = 5, align = "v", #labels = "auto",
               labels = c("(a) Kelch 13", "(b) Pfcrt-K76T", "(c) Pfmdr1-N86Y",
                          "(d) Pfmdr1-Y184F", "(e) Pfmdr1-D1246Y"),
               label_fontface = "plain",
               label_y = 0.95, label_x = 0.07, hjust = 0) %>%
  plot_grid(legs, nrow = 1, rel_widths = c(1, 0.2)) +
  theme(plot.margin = margin(t = 0, r = 0, b = 20, l = 25))

legs <- get_legend(p1)

ggdraw(p) +
  draw_label(
    "Record index",
    x = 0.4, y = 0.02,
    hjust = 0.5, vjust = 0
  ) +
  draw_label(
    "Recorded molecular marker prevalence",
    x = 0.02, y = 0.5,
    angle = 90,
    hjust = 0.5, vjust = 1
  )
ggsave("figures/crints_bb.png", scale = 1, height = 10, width = 9)


png("~/Desktop/Test.png", width = 2000, height = 2000)
par(mfrow = c(2,2))
acf(0, ylab = "xxx")
acf(0)
acf(0)
acf(0)
dev.off()

##############################################################################
# a plot of all preds in all years

preds <- rast("output/k13_marcse/gneiting_sparse/preds_all.grd")
coords <- xyFromCell(preds, cells(preds))
vals <- terra::extract(preds, coords)
df <- cbind(coords, vals) %>%
  pivot_longer(starts_with("2"),
               names_to = "lyr",
               values_to = "val") %>%
  mutate(year = substr(lyr, 1, 4),
         tag = substr(lyr, 11, 14),
         marker = str_extract(lyr, "[^_]+$")) %>%
  filter(tag == "medi") # pick out year and thingo

ggplot() +
  geom_sf(data = st_as_sf(afr), fill = "white") +
  geom_tile(aes(x = x, y = y, fill = val), data = df) +
  scale_fill_viridis_c(na.value = NA, "Prevalence", trans = "sqrt") +
  facet_wrap(~year, ncol = 6) +
  scale_x_continuous(breaks = seq(-20, 40, 20), "Longitude") +
  scale_y_continuous(breaks = seq(-20, 40, 20), "Latitude") +
  xlab("Longitude") +
  ylab("Latitude")

ggsave("figures/k13_all_years.png", scale = 1.7, height = 4, width = 5)

preds <- rast("output/k13_marcse/bb_gne/preds_all.tif")
coords <- xyFromCell(preds, cells(preds))
vals <- terra::extract(preds, coords)
df <- cbind(coords, vals) %>%
  pivot_longer(starts_with("2"),
               names_to = "lyr",
               values_to = "val") %>%
  mutate(year = substr(lyr, 1, 4),
         tag = substr(lyr, 11, 14),
         marker = str_extract(lyr, "[^_]+$")) %>%
  filter(marker == "50") # pick out year and thingo

ggplot() +
  geom_sf(data = st_as_sf(afr), fill = "white") +
  geom_tile(aes(x = x, y = y, fill = val), data = df) +
  scale_fill_viridis_c(na.value = NA, "Prevalence", trans = "sqrt") +
  facet_wrap(~year, ncol = 6) +
  scale_x_continuous(breaks = seq(-20, 40, 20), "Longitude") +
  scale_y_continuous(breaks = seq(-20, 40, 20), "Latitude") +
  xlab("Longitude") +
  ylab("Latitude")

ggsave("figures/k13_all_years_bb.png", scale = 1.7, height = 4, width = 5)

################################################################################
# let's crack into some partial dependence plots

marker <- "k13_marcse"
plots <- list(pfpr = plot_pdps(out_dir = paste0("output/", marker, "/bb_gne/"), 
                                         target = "pfpr"),
                        year = plot_pdps(out_dir = paste0("output/", marker, "/bb_gne/"), 
                                         target = "pfpr"))

plots <- read_rds("output/k13_marcse/pdps_plotted.rds")

for

################################################################################

cat(reports, file = to_report, append = TRUE, sep = "\n")







