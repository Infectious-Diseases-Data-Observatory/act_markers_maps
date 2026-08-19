# merge IDDO Surveyor data in data/raw/ with MARCSE dashboard data 
library(tidyverse)
library(readxl)
library(viridisLite)

# this script does all of our packages and brings in the WHO markers list:
source("code/data_clean/format_moldm.R")
to_report <- "output/stats_to_report.txt"
reports <- "#### K13 data read-in ####"

moldm <- format_moldm_k13("data/raw/db_20260211/novartis.csv", 
                          report_path = to_report) %>%
  mutate(Marker = strip_marker) %>%
  clean_up_pmids() %>%
  mutate(Year.Published = as.numeric(Year.Published)) %>%
  suppressWarnings()

# this is sitting in a separate repo at lu-harr/MARC_SEA_dashboard/
# put in a pull request to Stephanie-van-Wyk
marcse <- read.csv("../MARC_SEA_dashboard/tidied_k13_dashboard_data.csv") %>%
  # this is initialised in the format_moldm script:
  left_join(marker_reference, 
            by = join_by(Marker == marker))

# let's check for weirdness here:
marcse %>%
  filter(!PubMedID %in% moldm$PubMedID) %>%
  dplyr::select(PubMedID) %>%
  unique() %>%
  as.vector()
marcse %>%
  filter(str_length(PubMedID) != 8) %>%
  dplyr::select(PubMedID) %>%
  unique() %>%
  as.vector()

# check for weirdness here:  
moldm %>%
  filter(!PubMedID %in% marcse$PubMedID) %>%
  dplyr::select(PubMedID) %>%
  unique() %>%
  as.vector()

# have paid attention to PMIDs as joining variable in the various tidy-ups ...
to_add <- anti_join(marcse, moldm, by = join_by(PubMedID)) %>%
  filter(PubMedID != "Already in moldm") %>%
  # add these back in:
  bind_rows(filter(marcse, PubMedID == "Unpublished"))
reports <- c(reports,
             paste0("MARCSE Studies to add: ", length(unique(to_add$Title))))
reports <- c(reports,
             paste0("with Publication years: ",
             range(as.numeric(to_add$Year.Published), na.rm=TRUE),
             collapse = ","))

reports <- c(reports,
             paste0("with Patients: ", to_add %>%
                 group_by(Longitude, Latitude, year, Title, Tested) %>%
                 summarise(n = n()) %>%
                 ungroup() %>%
                 dplyr::select(Tested) %>%
                 sum() %>%
                 suppressMessages()))
reports <- c(reports,
             paste0("Or more conservatively: ", to_add %>%
                 group_by(Longitude, Latitude, year, Title) %>%
                 summarise(n = length(unique(Tested)), Tested = max(Tested)) %>%
                 ungroup() %>%
                 dplyr::select(Tested) %>%
                 sum() %>%
                 suppressMessages()))

moldm <- bind_rows(moldm,
                   anti_join(marcse, moldm, by = join_by(PubMedID)) %>%
                   bind_rows(filter(marcse, PubMedID == "Unpublished"))
                   ) %>%
  filter(PubMedID != "Already in moldm") %>% # a sneaky preprint snuck through
  mutate(mutant = !is.na(status))

# have a look for 535K ...
# moldm %>% 
#   filter(grepl("535", Marker)) %>%
#   dplyr::select(Country, District, Site.Name, Longitude, Latitude, Start.Year, End.Year,
#                 Marker, Present, Tested, PubMedID, Year.Published, Title, Authors, Publication.URL, from) %>%
#   write.csv("../k13_parasite_clearance/535_for_Stephanie.csv", row.names = FALSE)

# clean these up beforehand if I can:
moldm %>% 
  filter(Tested < Present) %>% 
  dplyr::select(Title, PubMedID) %>% 
  unique()

message("From here down it's pretty much as in the format_moldm script")

reports <- c(reports,
             paste0("joint Number total: ", 
                    moldm %>%
                      group_by(Longitude, Latitude, year, PubMedID) %>% 
                      summarise(n = length(unique(Tested)), Tested = max(Tested), Present = sum(Present)) %>%
                      ungroup() %>%
                      dplyr::select(Present) %>%
                      sum() %>%
                      suppressMessages()))

reports <- c(reports,
             paste0("joint Number of wildtypes: ", 
                    moldm %>%
                      filter(Marker == "wildtype") %>%
                      group_by(Longitude, Latitude, year, PubMedID) %>% 
                      summarise(n = length(unique(Tested)), Tested = max(Tested), Present = sum(Present)) %>%
                      ungroup() %>%
                      dplyr::select(Present) %>%
                      sum() %>%
                      suppressMessages()))

reports <- c(reports,
             paste0("joint Number of validated mutants: ", 
                    moldm %>%
                      filter(status == "Validated marker (Epi., Lab. & Clin.)") %>%
                      group_by(Longitude, Latitude, year, PubMedID) %>% 
                      summarise(n = length(unique(Tested)), Tested = max(Tested), Present = sum(Present)) %>%
                      ungroup() %>%
                      dplyr::select(Present) %>%
                      sum() %>%
                      suppressMessages()))

reports <- c(reports,
             paste0("joint Number of candidate mutants: ", 
                    moldm %>%
                      filter(status == "Candidate marker (Epi. & Clin.)") %>%
                      group_by(Longitude, Latitude, year, PubMedID) %>% 
                      summarise(n = length(unique(Tested)), Tested = max(Tested), Present = sum(Present)) %>%
                      ungroup() %>%
                      dplyr::select(Present) %>%
                      sum() %>%
                      suppressMessages()))


write.csv(moldm %>%
            filter(Present / Tested <= 1 & Tested > MIN_SAMPLE_SIZE), 
          "data/clean/moldm_marcse_with_markers.csv")

mutants <- moldm %>%
  filter(mutant) %>%
  filter(Present > 0) %>%
  group_by(Longitude, Latitude, year, Tested, Site.Name, Country) %>%
  summarise(Present = sum(Present), 
            pubs = paste0(unique(PubMedID), collapse=",")) %>%
  mutate(Site.Name = gsub(",", "", Site.Name)) %>%
  arrange(year) %>%
  dplyr::select(Longitude, Latitude, year, Tested, Present, Site.Name, Country, pubs) %>%
  ungroup() %>%
  suppressMessages()

reports <- c(reports,
             paste("Number of mutants on the list in the dataset:", 
                   moldm %>% 
                     filter(mutant) %>% 
                     dplyr::select(Marker) %>% 
                     ungroup() %>%
                     unique() %>% 
                     unlist() %>% 
                     length()))

reports <- c(reports,
             "TF-associated mutations in dataset: ")
reports <- c(reports,
             moldm %>% 
  filter(mutant) %>% 
  group_by(Marker) %>%
  summarise(n_present = sum(Present), n_tested = sum(Tested)) %>%
  filter(n_present > 0) %>%
  full_join(marker_reference, join_by(Marker == marker)) %>%
  arrange(desc(n_present)) %>%
  filter(!is.na(n_present)) %>% 
  as.data.frame() %>%
    do.call(what = paste, .))

wildtypes <- moldm %>%
  filter(Marker == "wildtype") %>%
  group_by(Longitude, Latitude, year, Tested, Site.Name, Country) %>%
  summarise(Present = sum(Present), n=n(), 
            pubs = paste0(unique(PubMedID), collapse=",")) %>% 
  # check how many simultaneous wildtype entries we have?
  mutate(Site.Name = gsub(",", "", Site.Name)) %>%
  arrange(year) %>%
  dplyr::select(Longitude, Latitude, year, Tested, Site.Name, Country, Present, pubs) %>%
  ungroup() %>%
  suppressMessages()

# poss duplicates?
# wildtypes[which(wildtypes$Present/wildtypes$Tested > 1),
#           c("Longitude", "Latitude",  "year", "Tested", "Site.Name",
#             "Present", "pubs")]

wildtypes_to_add <- anti_join(
  # locations in `wildtypes` that do not occur in `mutants`,
  # paying attention to `Tested` but NOT to `pubs`
  wildtypes %>%
    dplyr::select(Longitude, Latitude, year, Tested, Site.Name, Country),
  mutants %>%
    dplyr::select(Longitude, Latitude, year, Tested, Site.Name, Country)) %>%
  mutate(Present = 0) %>%
  suppressMessages()

write.csv(wildtypes_to_add,
          "data/clean/moldm_marcse_wildtypes_to_add.csv",
          row.names = FALSE)

reports <- c(reports,
             paste("Number of rows of wildtypes to add:", nrow(wildtypes_to_add)))
# 536 - 413 == 123 added rows from marcse

with_wildtypes <- full_join(mutants, wildtypes_to_add) %>%
  filter(Tested > MIN_SAMPLE_SIZE) %>%
  suppressMessages()

# Make sure anything that pops up down here is cleaned up:
with_wildtypes %>%
  filter(Present / Tested > 1)
marcse %>%
  filter(Present / Tested > 1) %>%
  dplyr::select(PubMedID, Country, Site.Name, Present, Tested, Marker, Title)

write.csv(with_wildtypes %>%
            dplyr::select(-c(pubs)) %>%
            filter(Present / Tested <= 1),
          "data/clean/moldm_marcse_k13_nomarker.csv",
          row.names = FALSE)

# 119666.02 ??? caused by imputing testeds from present/prevalence ??
# (Although could have sworn )
# used to be 82041

###############################################################################
with_wildtypes <- read.csv("data/clean/moldm_marcse_k13_nomarker.csv")
wildtypes_to_add <- read.csv("data/clean/moldm_marcse_wildtypes_to_add.csv")
moldm <- read.csv("data/clean/moldm_marcse_with_markers.csv")

plot(with_wildtypes$year, with_wildtypes$Present/with_wildtypes$Tested, 
     xlab="Year", ylab="Prevalence")
with_wildtypes %>% filter(Present / Tested > 0.5)

denom <- moldm %>%
  group_by(Longitude, Latitude, year, PubMedID, Country) %>%
  summarise(Tested = max(Tested)) %>%
  ungroup() %>%
  dplyr::select(Tested) %>%
  sum()

nume <- moldm %>%
  group_by(Longitude, Latitude, year, PubMedID, Country) %>%
  summarise(Tested = max(Tested)) %>%
  ungroup() %>%
  filter(Country %in% c("Uganda")) %>%
  dplyr::select(Tested) %>%
  sum()

reports <- c(reports,
             paste0("Proportion of records in Uganda: ", nume/denom))

# tmp <- moldm %>%
#   mutate(Year.Published = as.numeric(Year.Published)) %>%
#   filter(Country == "South Africa"  & !is.na(status))
# ggplot()+
#   geom_sf(data = afr %>% filter(name == "South Africa")) +
#   geom_point(aes(x = Longitude, y = Latitude, col = Marker), data = tmp)

to_vis <- with_wildtypes %>% 
  mutate(year_bin = cut(year, 
                        breaks = c(min(year)-1, seq(2009, 2021, 3), max(year)))) %>%
  arrange(Present/Tested, Tested) %>%
  drop_na(Longitude, Latitude)

ggplot() + 
  geom_sf(data = afr, fill = "white") + 
  geom_point(data = filter(to_vis, Present == 0), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested, col = "grey50"),
             fill = "grey60",  pch = 21, alpha = 0.5, stroke = 0.2) +
  #new_scale_color() +
  geom_point(data = filter(to_vis, Present > 0), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested,
                           fill = Present / Tested),
             col = "grey50", pch=21, stroke = 0.2) +
  scale_color_manual(name = "Absence", values = c("grey30"), labels=c("")) +
  scale_fill_viridis_c(name = "Prevalence", trans = "sqrt") +
  scale_size_continuous(name = "Tested", range = c(0.2, 4), trans = "sqrt") +
  facet_wrap(~ year_bin, ncol=3) +
  labs(title = "Prevalence of k13 markers in Africa") +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_x_continuous(breaks = seq(-20, 40, 20)) +
  scale_y_continuous(breaks = seq(-20, 40, 20)) +
  theme_grey() 
# ggsave("figures/archive/abs_grey.png", height=3, width=5, scale=2)
# spooky warning here ...

markers <- moldm %>%
  filter(mutant) %>% # no filter on Tested
  group_by(year, Marker) %>%
  summarise(present = sum(Present)) %>%
  rename(marker = Marker)

# put test intensity in the background
bg <- moldm %>%
  group_by(Longitude, Latitude, year, PubMedID) %>%
  summarise(Tested = max(Tested)) %>%
  ungroup() %>%
  group_by(year) %>%
  summarise(Tested = sum(Tested))

# markers to show in top panel
# markers_panel_a <- markers %>%
#   ungroup() %>%
#   group_by(marker) %>%
#   summarise(maxp = max(present), sump = sum(present)) %>%
#   filter(maxp > 100) %>%
#   arrange(desc(sump)) %>%
#   ungroup() %>%
#   # need this for legend ordering ...
#   mutate(marker = factor(marker, levels = marker)) %>%
#   dplyr::select(marker) %>%
#   unique()

# markers to show in bottom set of panels
markers_panel_b <- markers %>% 
  ungroup() %>%
  group_by(marker) %>%
  summarise(n = sum(present)) %>%
  arrange(desc(n)) %>%
  dplyr::slice(1:5) %>%
  bind_rows(data.frame(marker = "Others", n=1))

npal = 6
# sty <- markers_panel_b %>%
#   mutate(col = factor(rep(1:npal, 5)[1:length(marker)]),
#          linet = factor(rep(c(1:2), each = npal)[1:length(marker)]))

df <- markers %>%
  ungroup() %>%
  # filter(marker %in% markers_panel_a$marker) %>%
  # mutate(marker = factor(marker, levels = markers_panel_a$marker)) %>%
  mutate(marker = ifelse(marker %in% markers_panel_b$marker,
                         marker,
                         "Others"),
         marker = factor(marker, levels = markers_panel_b$marker)) %>%
  summarise(.by = c(year, marker), present = sum(present)) %>%
  # left_join(sty, join_by(marker==marker)) %>%
  # mutate(group = interaction(col, linet, sep = " / ")) %>%
  filter(year > 2006)

# chop off <2010? fix colour palette some more?

# ggplot(markers %>%
#          filter(Marker %in% markers_keep$Marker)) +
#   geom_line(aes(x = year, y = present, colour = Marker)) +
#   labs(title="Most prevalent markers by year")

wts <- wildtypes_to_add %>%
  dplyr::select(Longitude, Latitude, year, Tested) %>%
  merge(markers_panel_b, by = NULL) %>%
  #rename(y = marker) %>% # what was I doing with this again?
  rename(Marker = marker) %>%
  mutate(Present = 0)

year_breaks <- c(YEAR_LOWER_BOUND, seq(2012, 2021, 3), 2024)
markers_disagg <- moldm %>%
  filter(mutant) %>%
  mutate(Marker = ifelse(Marker %in% markers_panel_b$marker,
                         Marker, "Others")) %>%
  filter(!(Marker == "Others" & Present == 0) & # don't want zeroes for all other markers, just WTs as BG
           Tested > MIN_SAMPLE_SIZE) %>% # don't want prevalence == 1, tested < 5 points
  filter(year > YEAR_LOWER_BOUND) %>%
  group_by(Marker, Longitude, Latitude, year, PubMedID) %>%
  summarise(Present = sum(Present), Tested = max(unique(Tested))) %>%
  bind_rows(wts) %>%
  mutate(year_bin = cut(year, 
                        breaks = year_breaks,
                        labels = paste(head(year_breaks, -1) + 1,
                                       tail(year_breaks, -1),
                                       sep = " - "))) %>%
  filter(Longitude > -30) %>% # remove ocean point that's a bit far away
  mutate(Marker = factor(Marker, levels = markers_panel_b$marker))

iddoblue <- iddoPal::iddo_palettes$iddo[1]

# Combine color and linetype into aesthetics
# I would now like to control the order in which the lines are plotted
# It shouldn't be that hard
# I shouldn't have to spend longer than 5 mins on stackoverflow
# It doesn't matter that much
# This is how we end up with crap graphs !!!!!!
# df <- df %>%
#   filter(year >= 2005) %>%
#   left_join(markers_panel_b, by = join_by(marker==marker)) %>%
#   arrange(desc(n))

bg <- bg %>%
  filter(year >= 2005)

bg_scale <- 25
bg_col <- "grey65"

message("Caution: setting bg_scale and y limits manually")
# p1 <- 
#   ggplot() +
#   geom_bar(data = bg, 
#            aes(x = year, y = Tested / bg_scale), 
#            stat = "identity", fill = "grey75") +
#   geom_line(data = df, linewidth = 0.7,
#             aes(x = year, y = present, group = marker, color = marker, 
#                 linetype = marker)) +
#   geom_point(data = df, 
#              aes(x = year, y = present, group = marker, color = marker)) +
#   labs(
#     title = "(a) Detected mutations by year",
#     color = "Marker",
#     linetype = "Marker"
#   ) +
#   scale_color_manual(values = rep(c(iddoblue, "#c7047c", viridis(4), "#E37210"), 2)) +
#   scale_linetype_manual(values = rep(1:2, each = 7)) +
#   scale_y_continuous(sec.axis = sec_axis(~.*bg_scale, name="Number of tests"),
#                      limits = c(0, 515)) +
#   theme_minimal() +
#   xlab("Year") +
#   ylab("Mutations detected") +
#   theme_bw() +
#   theme(axis.text = element_text(size = 10),
#         legend.text = element_text(size = 10),
#         legend.key.spacing.y = unit(-0.3, "lines"),
#         legend.box.margin = margin(50, 6, 6, 6),
#         axis.text.y.right = element_text(color = bg_col),
#         axis.title.y.right = element_text(color = bg_col),
#         axis.ticks.y.right = element_line(color = bg_col)) +
#   #legend.justification.right = "bottom") +
#   scale_x_continuous(breaks = 2008:2025, limits = c(2008, 2025)) 
# p1

p1 <- 
  ggplot() +
  geom_bar(data = bg, 
           aes(x = year, y = Tested / bg_scale), 
           stat = "identity", fill = "grey75") +
  # geom_line(data = df, linewidth = 0.7,
  #           aes(x = year, y = present, group = marker, color = marker, 
  #               linetype = marker)) +
  geom_bar(data = df, aes(x = year, y = present, fill = marker), 
           stat = "identity", width = 0.6) +
  labs(
    title = "(a) Detected mutations by year",
    fill = "Marker",
  ) +
  scale_fill_manual(values = rep(c(iddoblue, "#c7047c", viridis(4), "#E37210"), 2)) +
  # scale_linetype_manual(values = rep(1:2, each = 7)) +
  scale_y_continuous(sec.axis = sec_axis(~.*bg_scale, name="Number of tests"),
                     limits = c(0, 1200), breaks = seq(0, 1200, 200)) +
  theme_minimal() +
  xlab("Year") +
  ylab("Mutations detected") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.spacing.y = unit(-0.3, "lines"),
        legend.box.margin = margin(50, 6, 6, 6),
        axis.text.y.right = element_text(color = bg_col),
        axis.title.y.right = element_text(color = bg_col),
        axis.ticks.y.right = element_line(color = bg_col)) +
  #legend.justification.right = "bottom") +
  scale_x_continuous(breaks = 2008:2024, limits = c(2007.5, 2024.5))
p1

# what happened to 2008 testing?

p2 <- ggplot() + 
  geom_sf(data = afr, fill = "white") + 
  geom_point(data = filter(markers_disagg, Present == 0),
             mapping = aes(x = Longitude, y = Latitude,
                           size = Tested, col = "grey50"),
             fill = "grey70",  pch = 21, alpha = 0.7, stroke = 0.2) +
  geom_point(data = markers_disagg %>% 
               filter(Present > 0) %>%
               arrange(Present / Tested), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested,
                           fill = Present / Tested),
             col = "grey50", pch=21, stroke = 0.2, alpha = 0.7) +
  scale_color_manual(name = "", values = c("grey30"), labels=c("Absence"), guide = "none") +
  scale_fill_viridis_c(name = "Prevalence", trans = "sqrt", breaks = c(0.01, 0.05, 0.25, 0.5)) +
  scale_size_continuous(name = "Sample size", range = c(0.2, 6), trans = "sqrt",
                        breaks = c(10, 100, 1000, 3000)) +
  # allows labelling of rows and columns:
  facet_grid(rows=vars(year_bin), cols=vars(Marker)) +
  #facet_wrap(~ year_bin + Marker, drop=FALSE) +
  labs(title = "(b) Prevalence of k13 markers in Africa") +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_x_continuous(breaks = seq(-20, 40, 20)) +
  scale_y_continuous(breaks = seq(-20, 40, 20)) +
  theme_bw() +
  theme(plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "bottom",
        legend.title = element_text(hjust = 0.5),
        legend.spacing.x = unit(15, "lines"),
        panel.spacing = unit(0, "lines")) +
  guides(fill = guide_colourbar(title.position = "top",
                                theme = theme(legend.key.width = unit(10, "lines"))),
         size = guide_legend(title.position = "top"))

p2

# p2 +
#    geom_rect(data = data.frame(xmin = 60, xmax = 80, ymin = 0, ymax = 40, year_bin = "(2015,2018]",
#                                Marker = "R622I"), 
#              aes(xmin=xmin, xmax=xmax, ymin=ymin,  ymax=ymax), 
#              colour="grey10", fill="grey85", linewidth=0.3)

library(grid)
library(cowplot)
# rect <- rectGrob(
#   x = 0.922,
#   y = 0.27,
#   width = unit(0.2, "in"),
#   height = unit(0.2, "in"),
#   gp = gpar(fill = "grey70", colour="grey50", alpha = 0.5)
# )

rect <- rectGrob(
  # x = 0.915,
  # y = 0.27,
  x = 0.4, 0.05,
  width = unit(0.2, "in"),
  height = unit(0.2, "in"),
  gp = gpar(fill = "grey70", colour="grey50", alpha = 0.5)
)

p2 <- ggdraw(p2) +
  draw_grob(rect) + 
  # draw_label("Absence", x = 0.96, y = 0.27, size = 10)
  draw_label("Absence", x = 0.45, y = 0.05, size = 10)
#theme(#strip.background = element_blank(),
#strip.text.x = element_blank(),
#strip.placement = "outside") +
#ggsave("figs/markers_disagg.png", height = 4, width = 5, scale = 2)

plot_grid(p1, p2, ncol = 1, rel_heights = c(0.6, 2))

ggsave("figures/markers_disagg_marcse.png", height = 6, width = 5, scale = 2)

# prevalence range
reports <- c(reports,
             paste0("Prevalence range: ",
             markers_disagg %>% 
  filter(Present > 0) %>% 
  ungroup() %>%
  mutate(prev = Present/Tested) %>% 
  dplyr::select(prev) %>% 
  unlist() %>%
  range() %>%
  suppressMessages()))

cat(reports, file = to_report, append = TRUE, sep = "\n")

################################################################################
moldm <- read.csv("data/clean/moldm_marcse_with_markers.csv")
moldm %>% 
  group_by(Marker) %>%
  summarise(n = sum(Present)) %>%
  arrange(desc(n))

tmp <- moldm %>%
  filter(Marker == "R561H" & 
           Present > 0 &
           year > 2017)

ggplot() +
  geom_sf(data = afr) +
  geom_jitter(aes(x = Longitude, y = Latitude, col = PubMedID, size = Present),  
             size = 2, 
             alpha = 0.5,
             data = tmp)

tmp %>%
  group_by(PubMedID, Title, year, Authors, Country) %>%
  summarise(n = n()) %>%
  as.data.frame()

tmp %>%
  filter(PubMedID %in% c("40744004", "40744006"))

################################################################################

ggplot() +
  geom_bar(data = bg, aes(x = year, y = Tested / bg_scale), 
           stat = "identity", fill = "grey75") +
  geom_line(data = df, size = 0.7,
            aes(x = year, y = present, group = marker, color = marker, linetype = marker)) +
  geom_point(data = df, 
             aes(x = year, y = present, group = marker, color = marker)) +
  labs(
    title = "Detected mutations by year - IDDO Surveyor + MARCSE-Africa data",
    color = "Marker",
    linetype = "Marker"
  ) +
  scale_color_manual(values = rep(c(viridis(4), "#E37210", iddoblue, "#c7047c"), 2)) +
  scale_linetype_manual(values = rep(1:2, each = 7)) +
  scale_y_continuous(sec.axis = sec_axis(~.*bg_scale, name="Number of tests"),
                     limits = c(0, 515)) +
  theme_minimal() +
  xlab("Year") +
  ylab("Mutations detected") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.spacing.y = unit(-0.3, "lines"),
        legend.box.margin = margin(50, 6, 6, 6),
        axis.text.y.right = element_text(color = bg_col),
        axis.title.y.right = element_text(color = bg_col),
        axis.ticks.y.right = element_line(color = bg_col)) +
  #legend.justification.right = "bottom") +
  scale_x_continuous(breaks = seq(2005, 2025, 2))
ggsave("~/Desktop/presentations/MARCSE/moldm_marcse_k13_time.png", scale = 0.8, width = 10, height = 5)



dat <- read.csv("data/clean/moldm_marcse_k13_nomarker.csv") %>%
  filter(year >= 2000) %>%
  mutate(year_bin = cut(year, breaks = seq(2000, 2025, 5))) %>%
  arrange(Present/Tested)

ggplot() + 
  geom_sf(data = afr, fill = "white") + 
  geom_point(data = filter(dat, Present == 0), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested, col = "grey50"),
             fill = "grey60",  pch = 21, alpha = 0.5, stroke = 0.2) +
  geom_point(data = filter(dat, Present > 0), 
             mapping = aes(x = Longitude, y = Latitude, 
                           size = Tested,
                           fill = Present / Tested),
             col = "grey50", pch=21, stroke = 0.2) +
  scale_color_manual(name = "Absence", values = c("grey30"), labels=c("")) +
  scale_fill_viridis_c(name = "Prevalence", trans = "sqrt") +
  scale_size_continuous(name = "Tested", range = c(0.5, 6), trans = "sqrt") +
  #facet_wrap(~ year_bin, ncol=3) +
  labs(title = "Prevalence of Kelch 13 markers - with unpublished records") +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_x_continuous(breaks = seq(-20, 40, 20)) +
  scale_y_continuous(breaks = seq(-20, 40, 20)) +
  theme_grey() 
ggsave("~/Desktop/presentations/MARCSE/moldm_marcse_k13.png", height = 6, width = 6)

# moldm %>%
#   filter(Marker == "P441L" & Present > 0) %>%
#   dplyr::select(Country, year, Tested, Present, Prevalence, PubMedID, from, Site.Name) %>%
#   mutate(Prevalence = Present / Tested) %>%
#   arrange(year, Prevalence)

# looks like some of the presences and tested are around the wrong way?
# Eloff 2025 entered in moldm twice :/


################################################################################
# 
# afr_with_centroids <- afr %>%
#   st_centroid() %>%
#   st_coordinates() %>%
#   bind_cols(afr) %>%
#   mutate(region = case_when(Y < -10 ~ "south",
#                             X > 22 & Y > -10 & Y < 22 ~ "east",
#                             X < 22 & Y > -10 & Y < 22 ~ "west")) %>%
#   st_as_sf()
# 
# safr <- afr_with_centroids %>%
#   filter(Y < -10) %>%
#   st_as_sf()
# 
# eafr <- afr_with_centroids %>%
#   filter(X > 22 & Y > -10 & Y < 22) %>%
#   st_as_sf()
# 
# wafr <- afr_with_centroids %>%
#   filter(X < 22 & Y > -10 & Y < 22) %>%
#   st_as_sf()
# 
# dat_tagged <- markers_disagg %>%
#   st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(afr))
# 
# dat_tagged <- dat_tagged %>%
#   mutate(eafr = st_intersects(., eafr) %>% lengths(),
#          wafr = st_intersects(., wafr) %>% lengths(),
#          safr = st_intersects(., safr) %>% lengths(),
#          region = case_when(eafr == 1 ~ "east",
#                             wafr == 1 ~ "west",
#                             safr == 1 ~ "south",
#                             TRUE ~ NA)) %>%
#   mutate(year_bin = cut(year, breaks = c(min(year) - 1, 2015, 2018, 2021, 2024)))
# 
# p <- ggplot() + 
#   geom_sf(data = afr_with_centroids, fill = "white") + 
#   geom_sf(data = filter(dat_tagged, Present == 0),
#           mapping = aes(size = Tested, col = "grey50"),
#           fill = "grey70",  pch = 21, alpha = 0.2, stroke = 0.2) +
#   geom_sf(data = dat_tagged %>% 
#             filter(Present > 0) %>%
#             arrange(Present / Tested), 
#           mapping = aes(size = Tested, fill = Present / Tested),
#           col = "grey50", pch=21, stroke = 0.2) +
#   scale_color_manual(name = "", values = c("grey30"), 
#                      labels=c("Absence"), guide = "none") +
#   scale_fill_viridis_c(name = "Prevalence", 
#                        trans = "sqrt",
#                        limits = c(min(dat_tagged$Present / dat_tagged$Tested),
#                                   max(dat_tagged$Present / dat_tagged$Tested))) +
#   scale_size_continuous(name = "Sample size", 
#                         range = c(0.2, 5), 
#                         trans = "sqrt") +
#   # allows labelling of rows and columns:
#   facet_grid(region ~ year_bin, scales = "free") +
#   #labs(title = "(b) Prevalence of k13 markers in Africa") +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   scale_x_continuous(breaks = seq(-20, 40, 20)) +
#   scale_y_continuous(breaks = seq(-20, 40, 20)) +
#   theme(plot.background = element_rect(fill='transparent', color=NA),
#         strip.text = element_blank())
# 
# p
# 
# p_safr <- ggplot() + 
#   geom_sf(data = safr, fill = "white") + 
#   geom_sf(data = filter(dat_tagged, Present == 0 & safr == 1),
#              mapping = aes(size = Tested, col = "grey50"),
#              fill = "grey70",  pch = 21, alpha = 0.2, stroke = 0.2) +
#   geom_sf(data = dat_tagged %>% 
#                filter(Present > 0 & safr == 1) %>%
#                arrange(Present / Tested), 
#              mapping = aes(size = Tested, fill = Present / Tested),
#              col = "grey50", pch=21, stroke = 0.2) +
#   scale_color_manual(name = "", values = c("grey30"), 
#                      labels=c("Absence"), guide = "none") +
#   scale_fill_viridis_c(name = "Prevalence", 
#                        trans = "sqrt",
#                        limits = c(min(dat_tagged$Present / dat_tagged$Tested),
#                                   max(dat_tagged$Present / dat_tagged$Tested))) +
#   scale_size_continuous(name = "Sample size", 
#                         range = c(0.2, 5), 
#                         trans = "sqrt") +
#   # allows labelling of rows and columns:
#   #facet_grid(rows=vars(year_bin), cols=vars(Marker)) +
#   #facet_wrap(~ year_bin + Marker, drop=FALSE) +
#   facet_wrap(~year_bin, ncol = 1) +
#   #labs(title = "(b) Prevalence of k13 markers in Africa") +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   scale_x_continuous(breaks = seq(-20, 40, 20)) +
#   scale_y_continuous(breaks = seq(-20, 40, 20)) +
#   theme(plot.background = element_rect(fill='transparent', color=NA),
#         strip.text = element_blank())
# 
# p_wafr <- ggplot() + 
#   geom_sf(data = wafr, fill = "white") + 
#   geom_sf(data = filter(dat_tagged, Present == 0 & wafr == 1),
#           mapping = aes(size = Tested, col = "grey50"),
#           fill = "grey70",  pch = 21, alpha = 0.2, stroke = 0.2) +
#   geom_sf(data = dat_tagged %>% 
#             filter(Present > 0 & wafr == 1) %>%
#             arrange(Present / Tested), 
#           mapping = aes(size = Tested, fill = Present / Tested),
#           col = "grey50", pch=21, stroke = 0.2) +
#   scale_color_manual(name = "", values = c("grey30"), 
#                      labels=c("Absence"), guide = "none") +
#   scale_fill_viridis_c(name = "Prevalence", 
#                        trans = "sqrt", 
#                        limits = c(min(dat_tagged$Present / dat_tagged$Tested),
#                                 max(dat_tagged$Present / dat_tagged$Tested))) +
#   scale_size_continuous(name = "Sample size", 
#                         range = c(0.2, 5), 
#                         trans = "sqrt") +
#   facet_wrap(~year_bin, ncol = 1) +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   scale_x_continuous(breaks = seq(-20, 40, 20)) +
#   scale_y_continuous(breaks = seq(-20, 40, 20)) +
#   theme(plot.background = element_rect(fill='transparent', color=NA),
#         strip.text = element_blank())
# 
# p_eafr <- ggplot() + 
#   geom_sf(data = eafr, fill = "white") + 
#   geom_sf(data = filter(dat_tagged, Present == 0 & eafr == 1),
#           mapping = aes(size = Tested, col = "grey50"),
#           fill = "grey70",  pch = 21, alpha = 0.2, stroke = 0.2) +
#   geom_sf(data = dat_tagged %>% 
#             filter(Present > 0 & eafr == 1) %>%
#             arrange(Present / Tested), 
#           mapping = aes(size = Tested, fill = Present / Tested),
#           col = "grey50", pch=21, stroke = 0.2) +
#   scale_color_manual(name = "", values = c("grey30"), 
#                      labels=c("Absence"), guide = "none") +
#   scale_fill_viridis_c(name = "Prevalence", 
#                        trans = "sqrt", 
#                        limits = c(min(dat_tagged$Present / dat_tagged$Tested),
#                                   max(dat_tagged$Present / dat_tagged$Tested))) +
#   scale_size_continuous(name = "Sample size", 
#                         range = c(0.2, 5), 
#                         trans = "sqrt") +
#   facet_wrap(~year_bin, ncol = 1) +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   scale_x_continuous(breaks = seq(-20, 40, 20)) +
#   scale_y_continuous(breaks = seq(-20, 40, 20)) +
#   theme(plot.background = element_rect(fill='transparent', color=NA),
#         strip.text = element_blank())
# p_eafr
# 
# plot_grid(p_wafr, p_eafr, p_safr, nrow = 1)
