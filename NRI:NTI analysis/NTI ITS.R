library(ggplot2)
library(tidyr)
library(dplyr)
library(ape)
library(phytools)
library(picante)

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

rhizo.shared = rhizo.norm[!(rhizo.norm$feature_id %in% rhizo.specialist$feature_id), ]
det.shared = det.norm[!(det.norm$feature_id %in% det.specialist$feature_id), ]
rd.shared = rd.norm[!(rd.norm$feature_id %in% rd.specialist$feature_id), ]

rhizo.specialist$type = "specialist"
det.specialist$type = "specialist"
rd.specialist$type = "specialist"

rhizo.shared$type = "shared"
det.shared$type = "shared"
rd.shared$type = "shared"

rhizo.types = rbind(rhizo.shared, rhizo.specialist)
det.types = rbind(det.shared, det.specialist)
rd.types = rbind(rd.shared, rd.specialist)
abund.data = read.csv("/Users/mfoley2/Documents/growth ecology paper/NRI NTI analysis/tube_level_counts_ITS.csv",header = TRUE)

abund.data = abund.data %>% select(-X)
abund.data = gather(abund.data, feature_id, abundance, 2:length(abund.data)) %>% rename(tube = sample)

levels(as.factor(abund.data$tube))

rhizo = abund.data %>% filter(tube == "193"| tube =="194"| tube =="195"|tube =="196"| tube =="197"|
                                tube =="198"|
                                tube =="199"|
                                tube =="200")

det = abund.data %>% filter(tube == "203"|
                              tube =="204"|
                              tube =="205"|
                              tube =="206"|
                              tube =="207"|
                              tube =="208")

rhizo.det = abund.data %>% filter(tube == "209"|
                                    tube =="210"|
                                    tube =="211"|
                                    tube =="212"|
                                    tube =="213"|
                                    tube =="214"|
                                    tube =="215"|
                                    tube =="216")

rhizo.2 = merge(rhizo, rhizo.types, by = "feature_id")
det.2 = merge(det, det.types, by = "feature_id")
rhizo.det.2 = merge(rhizo.det, rd.types, by = "feature_id")

length(unique(rd.types$feature_id))
length(unique(rhizo.det.2$feature_id))

levels(as.factor(rhizo.2$type))

###rhizo
rhizo.otu.shared = rhizo.2 %>% filter(type == "shared") %>% select(-type)
rhizo.otu.unique = rhizo.2 %>% filter(type == "specialist")  %>% select(-type)

rhizo.otu.shared = spread(rhizo.otu.shared, feature_id, abundance)
rhizo.otu.unique = spread(rhizo.otu.unique, feature_id, abundance)

library(vegan)

rhizo.otu.shared.2 = decostand(rhizo.otu.shared[,2:length(rhizo.otu.shared)], method = "total")
rhizo.otu.unique.2 = decostand(rhizo.otu.unique[,2:length(rhizo.otu.unique)], method = "total")

sum(rhizo.otu.unique.2[1,])
sum(rhizo.otu.unique.2[,1])

### detritus
det.otu.shared = det.2 %>% filter(type == "shared") %>% select(-type)
det.otu.unique = det.2 %>% filter(type == "specialist")  %>% select(-type)

det.otu.shared = spread(det.otu.shared, feature_id, abundance)
det.otu.unique = spread(det.otu.unique, feature_id, abundance)

library(vegan)

det.otu.shared.2 = decostand(det.otu.shared[,2:length(det.otu.shared)], method = "total")
det.otu.unique.2 = decostand(det.otu.unique[,2:length(det.otu.unique)], method = "total")

sum(det.otu.shared.2[1,])
sum(det.otu.unique.2[,1])

### rhizo + detritus
rhizo.det.otu.shared = rhizo.det.2 %>% filter(type == "shared") %>% select(-type)
rhizo.det.otu.unique = rhizo.det.2 %>% filter(type == "specialist")  %>% select(-type)

rhizo.det.otu.shared = spread(rhizo.det.otu.shared, feature_id, abundance)
rhizo.det.otu.unique = spread(rhizo.det.otu.unique, feature_id, abundance)

rhizo.det.otu.shared.2 = decostand(rhizo.det.otu.shared[,2:length(rhizo.det.otu.shared)], method = "total")
rhizo.det.otu.unique.2 = decostand(rhizo.det.otu.unique[,2:length(rhizo.det.otu.unique)], method = "total")

sum(rhizo.det.otu.unique.2[1,])
sum(rhizo.det.otu.unique.2[,1])

###NRI NTI analysis
tree<-read.newick("/Users/mfoley2/Library/CloudStorage/OneDrive-UniversitédeLausanne/github/NRI:NTI/DRIPITS_seqs_2x2.fa.aln.tre")

rownames(rhizo.otu.shared.2) = rhizo.otu.shared$tube
rownames(rhizo.otu.unique.2) = rhizo.otu.unique$tube

rownames(det.otu.shared.2) = det.otu.shared$tube
rownames(det.otu.unique.2) = det.otu.unique$tube

rownames(rhizo.det.otu.shared.2) = rhizo.det.otu.shared$tube
rownames(rhizo.det.otu.unique.2) = rhizo.det.otu.unique$tube

#rhizo shared
picCleanComm <- match.phylo.comm(phy = tree, comm = rhizo.otu.shared.2)$comm
picCleanTree <- match.phylo.comm(phy = tree, comm = rhizo.otu.shared.2)$phy

cophenDist <- cophenetic.phylo(picCleanTree)

rhizo.generalist.mntd.out <- ses.mntd(picCleanComm, 
                       cophenDist, 
                       null.model = "taxa.labels", 
                       abundance.weighted = FALSE, 
                       runs = 999)
head(ses.mntd.out)

#rhizo unique
picCleanComm <- match.phylo.comm(phy = tree, comm = rhizo.otu.unique.2)$comm
picCleanTree <- match.phylo.comm(phy = tree, comm = rhizo.otu.unique.2)$phy

cophenDist <- cophenetic.phylo(picCleanTree)

rhizo.specialist.mntd.out <- ses.mntd(picCleanComm, 
                       cophenDist, 
                       null.model = "taxa.labels", 
                       abundance.weighted = FALSE, 
                       runs = 999)
head(ses.mntd.out)

#det shared
picCleanComm <- match.phylo.comm(phy = tree, comm = det.otu.shared.2)$comm
picCleanTree <- match.phylo.comm(phy = tree, comm = det.otu.shared.2)$phy

cophenDist <- cophenetic.phylo(picCleanTree)

det.generalist..mntd.out <- ses.mntd(picCleanComm, 
                       cophenDist, 
                       null.model = "taxa.labels", 
                       abundance.weighted = FALSE, 
                       runs = 999)
head(ses.mntd.out)

#det unique
picCleanComm <- match.phylo.comm(phy = tree, comm = det.otu.unique.2)$comm
picCleanTree <- match.phylo.comm(phy = tree, comm = det.otu.unique.2)$phy

cophenDist <- cophenetic.phylo(picCleanTree)

det.specialist..mntd.out <- ses.mntd(picCleanComm, 
                       cophenDist, 
                       null.model = "taxa.labels", 
                       abundance.weighted = FALSE, 
                       runs = 999)
head(ses.mntd.out)

#rhizo det shared
picCleanComm <- match.phylo.comm(phy = tree, comm = rhizo.det.otu.shared.2)$comm
picCleanTree <- match.phylo.comm(phy = tree, comm = rhizo.det.otu.shared.2)$phy

cophenDist <- cophenetic.phylo(picCleanTree)

rd.generalist.mntd.out <- ses.mntd(picCleanComm, 
                       cophenDist, 
                       null.model = "taxa.labels", 
                       abundance.weighted = FALSE, 
                       runs = 999)
head(ses.mntd.out)

#rhizo det unique
picCleanComm <- match.phylo.comm(phy = tree, comm = rhizo.det.otu.unique.2)$comm
picCleanTree <- match.phylo.comm(phy = tree, comm = rhizo.det.otu.unique.2)$phy

cophenDist <- cophenetic.phylo(picCleanTree)

rd.specialist.mntd.out <- ses.mntd(picCleanComm, 
                       cophenDist, 
                       null.model = "taxa.labels", 
                       abundance.weighted = FALSE, 
                       runs = 999)
head(ses.mntd.out)

NRI <- as.matrix(-1 * ((ses.mntd.out[,2] - ses.mntd.out[,3]) / ses.mntd.out[,4]))

rownames(NRI) <- row.names(ses.mntd.out)
colnames(NRI) <- "NRI"

head(NRI)
