data.18O = read.csv("/Users/mfoley2/Library/CloudStorage/OneDrive-UniversitédeLausanne/ITS_qsip_results_18O_outliers_removed_with_taxonomy.csv", as.is = TRUE)

data.18O.summary = data.18O %>% group_by(c.trt, water.trt, feature_id, Phylum, Class, Order, Family, Genus) %>% 
  summarise(mean_rgr = mean(rgr.per.day, na.rm = TRUE)) %>% filter(!is.na(Phylum))

library(tidyr)
#normal
data.18O.effects = data.18O.summary %>% filter(water.trt == "Normal")
data.18O.effects = spread(data.18O.effects, c.trt, mean_rgr)
data.18O.effects$detritus.effect.rgr = data.18O.effects$`Rhizosphere + detritus` - data.18O.effects$Rhizosphere
data.18O.effects$rhizo.effect.rgr = data.18O.effects$`Rhizosphere + detritus` - data.18O.effects$Detritus

data.18O.effects = na.omit(data.18O.effects)

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
custom_palette <- c("#1f77b4", "#A6CEE3","#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
      
                     "#7f7f7f", "#bcbd22", "#17becf", "#a55194")


levels(as.factor(data.18O.effects$Phylum))
data.18O.effects$Phylum = gsub("p__Anthophyta", "Anthophyta", data.18O.effects$Phylum)
data.18O.effects$Phylum = gsub("p__Ascomycota", "Ascomycota", data.18O.effects$Phylum)
data.18O.effects$Phylum = gsub("p__Basidiomycota", "Basidiomycota", data.18O.effects$Phylum)
data.18O.effects$Phylum = gsub("p__Chytridiomycota", "Chytridiomycota", data.18O.effects$Phylum)
data.18O.effects$Phylum = gsub("p__Coniferophyta", "Coniferophyta", data.18O.effects$Phylum)
data.18O.effects$Phylum = gsub("p__Glomeromycota", "Glomeromycota", data.18O.effects$Phylum)
data.18O.effects$Phylum = gsub("p__Kickxellomycota", "Kickxellomycota", data.18O.effects$Phylum)
data.18O.effects$Phylum = gsub("p__Mucoromycota", "Mucoromycota", data.18O.effects$Phylum)
data.18O.effects$Phylum = gsub("p__Rozellomycota", "Rozellomycota", data.18O.effects$Phylum)
data.18O.effects$Phylum = gsub("unassigned", "Unassigned", data.18O.effects$Phylum)

p2 = ggplot(data.18O.effects, aes(shape = water.trt, color = Phylum)) + 
  xlab("rhizosphere effect") +
  scale_color_manual(values = custom_palette) +
  ylab("detritrusphere effect") +
  geom_abline(intercept = -0.06, slope = 1.30, linetype = "solid", color = "black") +
    geom_abline(intercept = -0.05, slope = 1.17, linetype = "dashed", color = "blue", size = 0.5) +
    
  geom_point(size = 5, pch = 21, aes(x=rhizo.effect.rgr, 
                                     y=detritus.effect.rgr, fill = Phylum), alpha = 0.75) +
  theme_bw() + 
  ylim(c(-2,2)) +
  xlim(c(-2,2)) + ggtitle("ITS") +
  theme(legend.position = "right",
        legend.text = element_text(size = 16),
        strip.background = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        plot.title = element_blank()) +
  scale_fill_manual(values = custom_palette) +
  annotate("text", x = 2, y = -1.75, label = "Observed: R2 = 0.20", hjust = 1, size = 16 / .pt) +
  annotate("text", x = 2, y = -2, label = "Permutation test: P=0.04", hjust = 1, size = 16 / .pt) +
  annotate("segment", x = -2, xend = -1.4, y = 1.85, yend = 1.85, linetype = "solid", color = "black", linewidth = 0.5) +
  annotate("text", x = -1.3, y = 1.85, label = "Observed fit", hjust = 0, size = 16 / .pt) +
  annotate("segment", x = -2, xend = -1.4, y = 1.55, yend = 1.55, linetype = "dashed", color = "blue") +
  annotate("text", x = -1.3, y = 1.55, label = "Permutation test", hjust = 0, size = 16 / .pt)

library(patchwork)
p1+p2

