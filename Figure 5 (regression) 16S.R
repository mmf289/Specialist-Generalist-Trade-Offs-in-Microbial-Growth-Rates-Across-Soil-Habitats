library(dplyr)
library(lmodel2)
library(ggplot2)
library(tidyr)

data.18O = read.csv("/Users/mfoley2/Library/CloudStorage/OneDrive-UniversitédeLausanne/16S_qsip_results_18O_outliers_removed_with_taxonomy.csv", as.is = TRUE)

#formatting
data.18O.summary = data.18O %>% group_by(c.trt, water.trt, feature_id, phylum, class, order, family, genus)  %>%
  mutate(phylum = if_else(phylum == "unclassified_Bacteria", "Unassigned", phylum)) %>%
  mutate(phylum = if_else(phylum == "candidate division WPS-1", "Candidate division WPS-1", phylum)) %>%
  filter(phylum != "unclassified_root") %>% 
  filter(phylum != "Cyanobacteria/Chloroplast") %>% 
  filter(phylum != "Plantae") %>%
  filter(water.trt == "Normal") %>%
  summarise(mean_rgr = mean(rgr.per.day, na.rm = TRUE))

#compute rhizosphere and detritusphere effects
#REffect = RD - D
#DEffect = RD - R
data.18O.effects = data.18O.summary
data.18O.effects = spread(data.18O.effects, c.trt, mean_rgr)
data.18O.effects$detritus.effect.rgr = data.18O.effects$`Rhizosphere + detritus` - data.18O.effects$Rhizosphere
data.18O.effects$rhizo.effect.rgr = data.18O.effects$`Rhizosphere + detritus` - data.18O.effects$Detritus

data.18O.effects = na.omit(data.18O.effects)

#regression
data.18O.effects = na.omit(data.18O.effects)
lmodel2(detritus.effect.rgr ~ rhizo.effect.rgr, data = data.18O.effects) #OLS is type I regression, MA is type II, SMA is type II. Report the two-tailed p-value.

#Permutation analysis 
fit_obs = lmodel2(detritus.effect.rgr ~ rhizo.effect.rgr, data = data.18O.effects)
regression_obs = fit_obs$regression.results

r2_obs = (fit_obs$r)^2
b_obs = regression_obs$Slope[regression_obs$Method == "SMA"]

set.seed(123)
nperm = 9999
b_perm = numeric(nperm)
r2_perm = numeric(nperm)
intercept_perm = numeric(nperm)

for (i in seq_len(nperm)) {
  
  #permute R, D, and RD columns
  det_perm = sample(data.18O.effects$Detritus)
  rhizo_perm = sample(data.18O.effects$Rhizosphere)
  rd_perm = sample(data.18O.effects$`Rhizosphere + detritus`)
  
  #compute RE and DE with permuted values
  RE_perm = rd_perm - det_perm
  DE_perm = rd_perm - rhizo_perm
  
  #fit regressions
  fit_perm = lmodel2(DE_perm ~ RE_perm)
  regression_perm = fit_perm$regression.results
  
  r2_perm[i] = (fit_perm$r)^2
  b_perm[i] = regression_perm$Slope[regression_perm$Method == "SMA"]
  intercept_perm[i] = regression_perm$Intercept[regression_perm$Method == "SMA"]
}

#summary stats
median(b_perm)
median(intercept_perm)

#hypothesis tests
(sum(b_perm >= b_obs) + 1) / (nperm + 1) #upper-tail P asks whether the slope is higher than the artifact (support for generalists)
(sum(b_perm <= b_obs) + 1) / (nperm + 1) #lower-tail P asks whether the slope is weaker than the artifact or negative, Since the artifact is built into the math, I think it follows that a slope weaker than the permuted distribution is evidence of a tradeoff.

#hypothesis tests
(sum(r2_perm >= r2_obs) + 1) / (nperm + 1) #upper-tail P asks whether the correlation is tighter than the artifact (support for biological coupling of responses)
(sum(r2_perm <= r2_obs) + 1) / (nperm + 1) #lower-tail P asks whether the correlation is weaker than the artifact

#sanity check model fit
ggplot(data.18O.effects, aes(x = rhizo.effect.rgr, y = detritus.effect.rgr)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_function(fun = purrr::as_mapper(~ .x*regression_obs[1,]$Slope + regression_obs[1,]$Intercept), colour = "red", n = 100) +  #OLS
  geom_function(fun = purrr::as_mapper(~ .x*regression_obs[2,]$Slope + regression_obs[2,]$Intercept), colour = "blue", n = 100) + #MA
  geom_function(fun = purrr::as_mapper(~ .x*regression_obs[3,]$Slope + regression_obs[3,]$Intercept), colour = "green", n = 100) #SMA

#figure
custom_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
                    "#999999", "#FDCDAC", "#C7E9C0", "#FFD92F", "#CD9600", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "pink", 
                    "purple", "#00BFC4")

p1 = ggplot(data.18O.effects, aes(shape = water.trt, color = phylum)) + 
  xlab("rhizosphere effect") +
  scale_color_manual(values = custom_palette) +
  ylab("detritrusphere effect") +
  geom_abline(intercept = 0.24, slope = 0.65, linetype = "solid", color = "black", size = 0.5) +
    geom_abline(intercept = 0.24, slope = 0.66, linetype = "dashed", color = "blue") +
  geom_point(size = 5, pch = 21, aes(x=rhizo.effect.rgr, 
                                     y=detritus.effect.rgr, fill = phylum), alpha = 0.75) +
  theme_bw() + 
  ylim(c(-2,2)) +
  xlim(c(-2,2)) + ggtitle("16S") +
  scale_fill_manual(values = custom_palette) +
  theme(legend.position = "right",
        legend.text = element_text(size = 16),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.title = element_blank(),
        plot.title = element_blank()) +
  annotate("text", x = 2, y = -1.75, label = "Observed: R2 = 0.28", hjust = 1, size = 16 / .pt) +
  annotate("text", x = 2, y = -2, label = "Permutation test: P=0.43", hjust = 1, size = 16 / .pt) +
  annotate("segment", x = -2, xend = -1.4, y = 1.85, yend = 1.85, linetype = "solid", color = "black", linewidth = 0.5) +
  annotate("text", x = -1.3, y = 1.85, label = "Observed fit", hjust = 0, size = 16 / .pt) +
  annotate("segment", x = -2, xend = -1.4, y = 1.55, yend = 1.55, linetype = "dashed", color = "blue") +
  annotate("text", x = -1.3, y = 1.55, label = "Permutation test", hjust = 0, size = 16 / .pt)
p1
