library(dplyr)
library(lsmeans)
library(ggplot2)
library(ggsci)
library(patchwork)
library(tidyr)

data.18O = read.csv("/Users/mfoley2/Library/CloudStorage/OneDrive-UniversitédeLausanne/16S_qsip_results_18O_outliers_removed_with_taxonomy.csv", as.is = TRUE)

data.18O.summary = data.18O %>% group_by(c.trt, water.trt, feature_id) %>% 
  summarise(mean_rgr = mean(rgr.per.day, na.rm = TRUE))

rhizo.norm = data.18O.summary %>% filter(c.trt == "Rhizosphere") %>% filter(water.trt == "Normal") %>% 
  filter(mean_rgr > 0) %>% ungroup() %>% select(feature_id) %>% unique()

det.norm = data.18O.summary %>% filter(c.trt == "Detritus") %>% filter(water.trt == "Normal") %>% 
  filter(mean_rgr > 0) %>% ungroup() %>% select(feature_id) %>% unique()

rd.norm = data.18O.summary %>% filter(c.trt == "Rhizosphere + detritus") %>% filter(water.trt == "Normal") %>% 
  filter(mean_rgr > 0) %>% ungroup() %>% select(feature_id) %>% unique()

rhizo.specialist = rhizo.norm[!(rhizo.norm$feature_id %in% det.norm$feature_id), ]
rhizo.specialist = rhizo.specialist[!(rhizo.specialist$feature_id %in% rd.norm$feature_id), ]

rhizo.specialist = rhizo.norm[!(rhizo.norm$feature_id %in% det.norm$feature_id), ]
rhizo.specialist = rhizo.specialist[!(rhizo.specialist$feature_id %in% rd.norm$feature_id), ]


det.specialist = det.norm[!(det.norm$feature_id %in% rhizo.norm$feature_id), ]
det.specialist = det.specialist[!(det.specialist$feature_id %in% rd.norm$feature_id), ]

rd.specialist = rd.norm[!(rd.norm$feature_id %in% rhizo.norm$feature_id), ]
rd.specialist = rd.specialist[!(rd.specialist$feature_id %in% det.norm$feature_id), ]

rhizo.generalist = rhizo.norm[!(rhizo.norm$feature_id %in% rhizo.specialist$feature_id), ]
det.generalist = det.norm[!(det.norm$feature_id %in% det.specialist$feature_id), ]
rd.generalist = rd.norm[!(rd.norm$feature_id %in% rd.specialist$feature_id), ]

rhizo.specialist$type = "specialist"
det.specialist$type = "specialist"
rd.specialist$type = "specialist"

rhizo.generalist$type = "generalist"
det.generalist$type = "generalist"
rd.generalist$type = "generalist"

rhizo.types = rbind(rhizo.generalist, rhizo.specialist)
det.types = rbind(det.generalist, det.specialist)
rd.types = rbind(rd.generalist, rd.specialist)

write.csv(rhizo.types, "/Users/mfoley2/Library/CloudStorage/OneDrive-UniversitédeLausanne/rhizosphere.microbes.csv")
write.csv(det.types, "/Users/mfoley2/Library/CloudStorage/OneDrive-UniversitédeLausanne/detritusphere.microbes.csv")
write.csv(rd.types, "/Users/mfoley2/Library/CloudStorage/OneDrive-UniversitédeLausanne/rhizosphere.detritusphere.microbes.types.csv")

rhizo.norm = data.18O %>% filter(c.trt == "Rhizosphere") %>% filter(water.trt == "Normal") %>% 
  filter(rgr.per.day > 0)

det.norm = data.18O %>% filter(c.trt == "Detritus") %>% filter(water.trt == "Normal") %>% 
  filter(rgr.per.day > 0)

rd.norm = data.18O %>% filter(c.trt == "Rhizosphere + detritus") %>% filter(water.trt == "Normal") %>% 
  filter(rgr.per.day > 0)

rhizo.norm = merge(rhizo.norm, rhizo.types, by = "feature_id")
det.norm = merge(det.norm, det.types, by = "feature_id")
rd.norm = merge(rd.norm, rd.types, by = "feature_id")

data = do.call(rbind, list(det.norm, rhizo.norm, rd.norm))
sum = data %>% group_by(water.trt, c.trt, type, tube) %>% summarise(mean_rgr = mean(rgr.per.day))

anova2<-aov(mean_rgr ~ c.trt*type, data= sum)
summary(anova2)

#post hoc test
sum.wide = spread(sum, type, mean_rgr)
ph.rhizo = sum.wide %>% filter(c.trt == "Rhizosphere")
ph.det = sum.wide %>% filter(c.trt == "Detritus")
ph.rd = sum.wide %>% filter(c.trt == "Rhizosphere + detritus")

t.test(ph.rhizo$specialist, ph.rhizo$generalist)
t.test(ph.det$specialist, ph.det$generalist)
t.test(ph.rd$specialist, ph.rd$generalist)

#compute SE
sum2 = sum %>% group_by(water.trt, c.trt, type) %>% summarise(mean = mean(mean_rgr), sd=sd(mean_rgr), n=n())
sum2$se = sum2$sd/sqrt(sum2$n)

sum2$c.trt = gsub("Detritus", "detritusphere", sum2$c.trt)
sum2$c.trt = gsub("Rhizosphere", "rhizosphere", sum2$c.trt)
sum2$c.trt = gsub("rhizosphere + detritus", "rhizosphere + detritusphere", sum2$c.trt,  fixed = T)

p1 = ggplot(sum2, aes(x = c.trt, y = mean, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("16S") +
  ylim(0,0.75) +
  labs(y=expression(paste("relative growth rate (day",""^{-1}, ")"))) +
  scale_fill_npg() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, mapping = TRUE), width = 0.2, position = position_dodge(width = 0.9)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 17),
        axis.text.y = element_text(size = 17),
        axis.title.x= element_blank(),
        axis.title.y=element_text(size = 17),
        panel.grid.minor=element_blank(),
        plot.title = element_text(hjust=0.5, size=18),
        panel.background = element_rect(color = "black", fill = "white"),
        strip.background = element_rect(color = "white", fill = "white"),
        strip.text = element_text(size = 17),
        legend.position = "none")
p1

##fungi
data.18O = read.csv("/Users/mfoley2/Library/CloudStorage/OneDrive-UniversitédeLausanne/ITS_qsip_results_18O_outliers_removed_with_taxonomy.csv", as.is = TRUE)

data.18O.summary = data.18O %>% group_by(c.trt, water.trt, feature_id) %>% 
  summarise(mean_rgr = mean(rgr.per.day, na.rm = TRUE))

rhizo.norm = data.18O.summary %>% filter(c.trt == "Rhizosphere") %>% filter(water.trt == "Normal") %>% 
  filter(mean_rgr > 0) %>% ungroup() %>% select(feature_id) %>% unique()

det.norm = data.18O.summary %>% filter(c.trt == "Detritus") %>% filter(water.trt == "Normal") %>% 
  filter(mean_rgr > 0) %>% ungroup() %>% select(feature_id) %>% unique()

rd.norm = data.18O.summary %>% filter(c.trt == "Rhizosphere + detritus") %>% filter(water.trt == "Normal") %>% 
  filter(mean_rgr > 0) %>% ungroup() %>% select(feature_id) %>% unique()

rhizo.specialist = rhizo.norm[!(rhizo.norm$feature_id %in% det.norm$feature_id), ]
rhizo.specialist = rhizo.specialist[!(rhizo.specialist$feature_id %in% rd.norm$feature_id), ]

rhizo.specialist = rhizo.norm[!(rhizo.norm$feature_id %in% det.norm$feature_id), ]
rhizo.specialist = rhizo.specialist[!(rhizo.specialist$feature_id %in% rd.norm$feature_id), ]

det.specialist = det.norm[!(det.norm$feature_id %in% rhizo.norm$feature_id), ]
det.specialist = det.specialist[!(det.specialist$feature_id %in% rd.norm$feature_id), ]

rd.specialist = rd.norm[!(rd.norm$feature_id %in% rhizo.norm$feature_id), ]
rd.specialist = rd.specialist[!(rd.specialist$feature_id %in% det.norm$feature_id), ]

rhizo.generalist = rhizo.norm[!(rhizo.norm$feature_id %in% rhizo.specialist$feature_id), ]
det.generalist = det.norm[!(det.norm$feature_id %in% det.specialist$feature_id), ]
rd.generalist = rd.norm[!(rd.norm$feature_id %in% rd.specialist$feature_id), ]

rhizo.specialist$type = "specialist"
det.specialist$type = "specialist"
rd.specialist$type = "specialist"

rhizo.generalist$type = "generalist"
det.generalist$type = "generalist"
rd.generalist$type = "generalist"

rhizo.types = rbind(rhizo.generalist, rhizo.specialist)
det.types = rbind(det.generalist, det.specialist)
rd.types = rbind(rd.generalist, rd.specialist)

write.csv(rhizo.types, "/Users/mfoley2/Library/CloudStorage/OneDrive-UniversitédeLausanne/rhizosphere.taxa.ITS.csv")
write.csv(det.types, "/Users/mfoley2/Library/CloudStorage/OneDrive-UniversitédeLausanne/detritusphere.taxa.ITS.csv")
write.csv(rd.types, "/Users/mfoley2/Library/CloudStorage/OneDrive-UniversitédeLausanne/rhizosphere.detritusphere.taxa.ITS.csv")

rhizo.norm = data.18O %>% filter(c.trt == "Rhizosphere") %>% filter(water.trt == "Normal") %>% 
  filter(rgr.per.day > 0)

det.norm = data.18O %>% filter(c.trt == "Detritus") %>% filter(water.trt == "Normal") %>% 
  filter(rgr.per.day > 0)

rd.norm = data.18O %>% filter(c.trt == "Rhizosphere + detritus") %>% filter(water.trt == "Normal") %>% 
  filter(rgr.per.day > 0)

rhizo.norm = merge(rhizo.norm, rhizo.types, by = "feature_id")
det.norm = merge(det.norm, det.types, by = "feature_id")
rd.norm = merge(rd.norm, rd.types, by = "feature_id")

data = do.call(rbind, list(det.norm, rhizo.norm, rd.norm))
sum = data %>% group_by(water.trt, c.trt, type, tube) %>% summarise(mean_rgr = mean(rgr.per.day))

anova2<-aov(mean_rgr ~ c.trt*type, data= sum)
summary(anova2)

#post hoc test
sum.wide = spread(sum, type, mean_rgr)
ph.rhizo = sum.wide %>% filter(c.trt == "Rhizosphere")
ph.det = sum.wide %>% filter(c.trt == "Detritus")
ph.rd = sum.wide %>% filter(c.trt == "Rhizosphere + detritus")

t.test(ph.rhizo$specialist, ph.rhizo$generalist)
t.test(ph.det$specialist, ph.det$generalist)
t.test(ph.rd$specialist, ph.rd$generalist)

sum2 = sum %>% group_by(water.trt, c.trt, type) %>% summarise(mean = mean(mean_rgr), sd=sd(mean_rgr), n=n())
sum2$se = sum2$sd/sqrt(sum2$n)

sum2$c.trt = gsub("Detritus", "detritusphere", sum2$c.trt)
sum2$c.trt = gsub("Rhizosphere", "rhizosphere", sum2$c.trt)
sum2$c.trt = gsub("rhizosphere + detritus", "rhizosphere + detritusphere", sum2$c.trt,  fixed = T)


p2 = ggplot(sum2, aes(x = c.trt, y = mean, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("ITS") +
  ylim(0,0.75) +
  labs(y=expression(paste("relative growth rate (day",""^{-1}, ")"))) +
  scale_fill_npg() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, mapping = TRUE), width = 0.2, position = position_dodge(width = 0.9)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 17),
        axis.text.y = element_blank(),
        axis.title.x= element_blank(),
        axis.title.y=element_blank(),
        panel.grid.minor=element_blank(),
        plot.title = element_text(hjust=0.5, size=18),
        panel.background = element_rect(color = "black", fill = "white"),
        strip.background = element_rect(color = "white", fill = "white"),
        strip.text = element_text(size = 17),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 17))
p2
p1 + p2

