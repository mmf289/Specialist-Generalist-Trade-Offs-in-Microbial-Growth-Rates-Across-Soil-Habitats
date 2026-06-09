library(dplyr)
library(ggtern)
library(plotly)
library(readr)
library(tidyr)
library(patchwork)
library(scales)

data = read.csv("/Users/mfoley2/Desktop/msystems resubmit/16S_qsip_results_18O_outliers_removed_with_taxonomy_strategy_annotated.csv", as.is = TRUE)

dat = data %>% group_by(c.trt, water.trt, feature_id) %>% 
  filter(phylum != "Plantae") %>% filter(phylum != "unclassified_root") %>%
  summarise(mean_rgr = mean(rgr.per.day, na.rm = TRUE)) %>% 
  filter(water.trt == "Normal") %>% ungroup() %>% select(feature_id, c.trt, mean_rgr)

dat$c.trt[dat$c.trt == "Detritus"] <- "Detritusphere"
dat$c.trt[dat$c.trt == "Rhizosphere + detritus"] <- "Rhizosphere + detritusphere"

dat_wide = spread(dat, c.trt, mean_rgr)
dat_wide[is.na(dat_wide)] <- 0

tax = data %>% select(phylum, feature_id)
tax$phylum <- gsub("unclassified_Bacteria", "unassigned", tax$phylum)

dat_wide_2 = merge(dat_wide, tax, by = "feature_id", all.x = T)

ggtern(data = dat_wide_2, aes(Detritusphere, Rhizosphere, `Rhizosphere + detritusphere`)) +
  theme_light() +
  theme_hidetitles() +
  theme_showarrows() +
  theme(tern.axis.text = element_text(size = 12),
        tern.axis.arrow = element_line(linewidth = 1, color = "black"),
        legend.position = "left",
        legend.text = element_text(size = 11)) +
  stat_hex_tern(
    bins  =15,                        # resolution of the tiling
    aes(fill = after_stat(count)),     # number of points in each hexagon
    alpha = 1, colour = NA
  )  + scale_fill_gradientn(colors=c("#ddf1da", "#abdda4","#e6f598", "#fee08b","#fdae61","#f46d43","#d53e4f"),
                            na.value="grey90",
                            limits  = c(0, 400),       # show 0–15; everything >15 is “15”
                            oob     = scales::squish,        # keep >15, just squash to top colour
                            name    = "ASV count")  

data.fun = read.csv("/Users/mfoley2/Documents/growth ecology paper/ITS_qsip_results_18O_outliers_removed_with_taxonomy_strategy_annotated.csv", as.is = TRUE)

###fungi
dat.fun = data.fun %>% group_by(c.trt, water.trt, feature_id) %>% 
  filter(!is.na(Phylum)) %>%
  summarise(mean_rgr = mean(rgr.per.day, na.rm = TRUE)) %>% 
  filter(water.trt == "Normal") %>% ungroup() %>% select(feature_id, c.trt, mean_rgr)

dat.fun$c.trt[dat.fun$c.trt == "Detritus"] <- "Detritusphere"
dat.fun$c.trt[dat.fun$c.trt == "Rhizosphere + detritus"] <- "Rhizosphere + detritusphere"


dat_wide_fun = spread(dat.fun, c.trt, mean_rgr)
dat_wide_fun[is.na(dat_wide_fun)] <- 0

tax = data.fun %>% select(Phylum, feature_id)

tax$Phylum <- gsub("p_", "", tax$Phylum)
tax$Phylum <- gsub("_", "", tax$Phylum)


dat_wide_fun_2 = merge(dat_wide_fun, tax, by = "feature_id", all.x = T)

ggtern(data=dat_wide_fun_2,
       aes(Detritusphere, Rhizosphere, `Rhizosphere + detritusphere`)) +
  theme_light() +
  theme_hidetitles() +
  theme_showarrows() +
  theme(tern.axis.text = element_text(size = 12),
        tern.axis.arrow = element_line(linewidth = 1, color = "black"),
        legend.position = "none") +
  stat_hex_tern(
    bins  =15,                        # resolution of the tiling
    aes(fill = after_stat(count)),     # number of points in each hexagon
    alpha = 1, colour = NA
  )  + scale_fill_gradientn(colors=c("#ddf1da", "#abdda4","#e6f598", "#fee08b"),
                            na.value="grey90",
                            limits  = c(0, 100),       # show 0–15; everything >15 is “15”
                            oob     = scales::squish,        # keep >15, just squash to top colour
                            name    = "ASV count")  

