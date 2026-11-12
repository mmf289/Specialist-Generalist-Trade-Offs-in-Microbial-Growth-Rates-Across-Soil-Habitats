devtools::install_github("jeffkimbrel/qSIP2")

library(dplyr)
library(ggplot2)
library(qSIP2)

packageVersion("qSIP2")

source_df<- read.csv("~/Desktop/H3_qsip_redo/jeff code/18O.lab.data.density.corrected.csv", as.is = TRUE)

#remove the following samples that were not sequenced so that Jeff's code runs smooth:
source_df = source_df %>% filter(sipID != 201) %>% filter(sipID != 202) %>% filter(sipID != 232) %>% #Michi is resequencing
  filter(sample.name != "195_13") %>% #not included in seq plans
  filter(sample.name != "195_14") %>% #not included in seq plans
  filter(sample.name != "195_15") %>% #not included in seq plans
  filter(sample.name != "195_16") %>% #not included in seq plans
  filter(sample.name != "195_17") %>% #not included in seq plans
  filter(sample.name != "195_13") %>% #not included in seq plans
  filter(sample.name != "197_10") %>% #lost during seq
  filter(sample.name != "220_17") %>% #not included in seq plans
  filter(sample.name != "220_18") %>% #not included in seq plans
  filter(sample.name != "228_16") %>% #lost during seq
  filter(sample.name != "235_10") #lost during seq


source_df_1 = source_df %>% select(sipID, c.trt, water.trt, harvest, isotope.trt) %>%
  unique()

source_object <- qsip_source_data(source_df_1,
                                  isotope = "isotope.trt",
                                  isotopolog = "c.trt",
                                  source_mat_id = "sipID"
)

sample_object <- qsip_sample_data(source_df,
                                  sample_id = "sample.name",
                                  source_mat_id = "sipID",
                                  gradient_position = "Fraction",
                                  gradient_pos_density = "density_corrected",
                                  gradient_pos_amt = "X16S.copies.uL"
)

asv_table<-read.delim("~/Desktop/H3_qsip/16S_data/DRIP16S_table_2x2.txt")

colnames(asv_table) <- gsub("X", "", colnames(asv_table))

asv_table_18O <- select(asv_table, -contains("F")) #remove 13C samples
asv_table_18O <- select(asv_table_18O , -contains("B")) #remove bulk samples
asv_table_18O <- select(asv_table_18O , -contains("U")) #remove unfractionated samples
asv_table_18O <- select(asv_table_18O , contains("_")) #remove additional unfractionated samples
asv_table_18O$ASV = asv_table$ASV

asv_table_18O <- select(asv_table_18O , -contains("195_19")) #outside of fraction range
asv_table_18O <- select(asv_table_18O , -contains("199_20")) #outside of fraction range
asv_table_18O <- select(asv_table_18O , -contains("204_20")) #outside of fraction range
asv_table_18O <- select(asv_table_18O , -contains("204_4")) #outside of fraction range
asv_table_18O <- select(asv_table_18O , -contains("204_5")) #outside of fraction range
asv_table_18O <- select(asv_table_18O , -contains("204_6")) #outside of fraction range
asv_table_18O <- select(asv_table_18O , -contains("204_1_")) #outside of fraction range
asv_table_18O <- select(asv_table_18O , -contains("204_2_")) #outside of fraction range
asv_table_18O <- select(asv_table_18O , -contains("204_211")) #outside of fraction range
asv_table_18O <- select(asv_table_18O , -contains("204_22")) #outside of fraction range
asv_table_18O <- select(asv_table_18O , -contains("204_3")) #outside of fraction range

feature_object <- qsip_feature_data(asv_table_18O,
                                    feature_id = "ASV"
)

qsip_object <- qsip_data(
  source_data = source_object,
  sample_data = sample_object,
  feature_data = feature_object
)

plot_source_wads(qsip_object, group = c("water.trt", "isotopolog")) 

plot_sample_curves(qsip_object)

show_comparison_groups(qsip_object,
                       group = c("water.trt", "isotopolog"),
                       isotope = "isotope",
                       source_mat_id = "source_mat_id"
)

###create comparisons
detritus_normal_204 <- run_feature_filter(qsip_object,
                                   unlabeled_source_mat_ids = c(203, 205, 207),
                                   labeled_source_mat_ids = 204,
                                   min_unlabeled_sources = 2,
                                   min_labeled_sources = 1,
                                   min_unlabeled_fractions = 2,
                                   min_labeled_fractions = 2
)


detritus_normal_204 <- run_resampling(detritus_normal_204,
                                      resamples = 1000,
                                      with_seed = 17,
                                      progress = FALSE,
                                      allow_failures = T,
)

plot_successful_resamples(detritus_normal_204)

detritus_normal_206 <- run_feature_filter(qsip_object,
                                          unlabeled_source_mat_ids = c(203, 205, 207),
                                          labeled_source_mat_ids = 206,
                                          min_unlabeled_sources = 2,
                                          min_labeled_sources = 1,
                                          min_unlabeled_fractions = 2,
                                          min_labeled_fractions = 2
)


detritus_normal_206 <- run_resampling(detritus_normal_206,
                                      resamples = 1000,
                                      with_seed = 17,
                                      progress = FALSE,
                                      allow_failures = T,
)

plot_successful_resamples(detritus_normal_206)

detritus_normal_208 <- run_feature_filter(qsip_object,
                                          unlabeled_source_mat_ids = c(203, 205, 207),
                                          labeled_source_mat_ids = 208,
                                          min_unlabeled_sources = 2,
                                          min_labeled_sources = 1,
                                          min_unlabeled_fractions = 2,
                                          min_labeled_fractions = 2
)


detritus_normal_208 <- run_resampling(detritus_normal_208,
                                      resamples = 1000,
                                      with_seed = 17,
                                      progress = FALSE,
                                      allow_failures = T,
)

plot_successful_resamples(detritus_normal_208)

detritus_normal_204 <- run_EAF_calculations(detritus_normal_204)
detritus_normal_206 <- run_EAF_calculations(detritus_normal_206)
detritus_normal_208 <- run_EAF_calculations(detritus_normal_208)

d_n_204 <- summarize_EAF_values(detritus_normal_204, confidence = 0.9) %>% mutate(tube = "204")
d_n_206 <- summarize_EAF_values(detritus_normal_206, confidence = 0.9) %>% mutate(tube = "206")
d_n_208 <- summarize_EAF_values(detritus_normal_208, confidence = 0.9) %>% mutate(tube = "208")

all.data = rbind(d_n_204, d_n_206)
all.data.2 = rbind(all.data, d_n_208)

ggplot(all.data.2) + 
  geom_boxplot(pch = 21, color = "black", aes(x=as.factor(tube), y=observed_EAF))

write.csv(all.data.2, "~/Desktop/H3_qsip_redo/jeff code/16S/detritus_normal.csv")

