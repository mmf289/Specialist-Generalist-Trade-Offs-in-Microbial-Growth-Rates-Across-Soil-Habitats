setwd("~/Desktop/H3_qsip_redo/jeff code")

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

data = read.csv("16S/results/qsip_results_18O.csv",header = TRUE)

#soil water enrichment was 0.693908785 atm % across all treatments
#assuming mu = 0.6, 0.693908785 * 0.6 = 0.416
data$eaf_max = 0.42

#visualize current eaf values - needs to be updated to show all 45 treatments
data$experimental_id = paste(data$c.trt, data$water.trt, sep="_")
sample_codes = levels(as.factor(data$experimental_id))

ggplot_list <- list()

for (name in sample_codes) {
  data.viz <- data %>% filter(experimental_id == name)
  p1 <- ggplot(data.viz, aes(x = observed_EAF)) +
    ggtitle(name) +
    geom_histogram() +
    geom_vline(xintercept = data.viz$eaf_max, color = "red") +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
          axis.title.x = element_text(size = 10),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(0.05, "inch"),
          legend.background = element_rect(fill = "transparent"))
  ggplot_list[[name]] <- p1
}

fig1 = ggplot_list[[1]] + ggplot_list[[2]] + ggplot_list[[3]] + ggplot_list[[4]] + ggplot_list[[5]] + ggplot_list[[6]]
fig1

####Assign taxa labels as too fast, negative, or normal
data.pos = data %>% filter(observed_EAF > 0) #for now, were only concerned with fast growing outliers so remove negatives
data.pos$difference = data.pos$eaf_max - data.pos$observed_EAF
too.fast = data.pos %>% filter(difference < 0) %>% select(-difference)
normal.eaf = data.pos %>% filter(difference > 0) %>% select(-difference)
negative = data %>% filter(observed_EAF<0)

too.fast$status = "too fast"
normal.eaf$status = "normal eaf"
negative$status = "negative"

#calculate sd
data.1 = rbind(too.fast, normal.eaf, negative)
data.1 = data.1 %>% group_by(feature_id, experimental_id) %>% mutate(sd = sd(observed_EAF, na.rm = T), n = n())

#visualize outliers on fast and slow end of spectrum
data.viz = data.1 %>% filter(status == "too fast")

ggplot(data.viz, aes(x=sd, y=observed_EAF)) + 
  geom_point(aes(color=status)) +facet_wrap(~experimental_id)

#sd of good data
data.1.b = data.1 %>% filter(status == "normal eaf") %>% group_by(feature_id, experimental_id) %>% mutate(sd = sd(observed_EAF), n = n())
ggplot(data.1.b, aes(y=sd, x = experimental_id)) + geom_boxplot()

#anything with SD above 0.15 among technical replicates will be removed
outliers = data.1 %>% filter(status == "negative" | status == "too fast") %>% filter(sd > 0.15) %>% mutate(remove.id = paste(feature_id, experimental_id, tube, sep="_"))
data.1$remove.id = paste(data.1$feature_id, data.1$experimental_id, data.1$tube, sep="_")
data_outliers_removed <- data.1[!data.1$remove.id %in% outliers$remove.id,]

#summary stats
length(unique(data_outliers_removed$feature_id))
data.summary = data_outliers_removed %>% group_by(experimental_id, feature_id) %>% summarise(mean = mean(observed_EAF))

#Visualize 
ggplot_list <- list()

# Loop through each value in sample_codes
for (name in sample_codes) {
  # Filter the data
  data.viz <- data_outliers_removed %>% filter(experimental_id == name)
  
  # Create the ggplot
  p1 <- ggplot(data.viz, aes(x = observed_EAF)) +
    ggtitle(name) +
    geom_histogram() +
    geom_vline(xintercept = data.viz$eaf_max, color = "red") +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
          axis.title.x = element_text(size = 10),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(0.05, "inch"),
          legend.background = element_rect(fill = "transparent"))
  
  # Store the ggplot in the list
  ggplot_list[[name]] <- p1
}

fig2 = ggplot_list[[1]] + ggplot_list[[2]] + ggplot_list[[3]] + ggplot_list[[4]] + ggplot_list[[5]] + ggplot_list[[6]] 

fig2

#correct negative values with bram's approach
prop_measures_for_correction = 0.10

#compute shifts for ember correction 
results_list <- list()

# Iterate through sample_names
for (name in sample_codes) {
  dat1 <- data_outliers_removed %>% filter(experimental_id == name)
  dat1 <- dat1[!is.na(dat1$observed_EAF),]
  dat1 <- dat1[order(dat1$observed_EAF),]
  dat2 <- dat1[1:floor(prop_measures_for_correction * nrow(dat1)),]
  shift2 <- median(dat2$observed_EAF)
  results_list[[name]] <- shift2
}

#merge results with original dataframe
my_output <- data.frame((sapply(results_list,c)))
my_output$experimental_id = rownames(my_output)
names(my_output) = c("shift", "experimental_id")

dat = merge(data_outliers_removed, my_output, by="experimental_id", all.x = TRUE)

#correct eaf values
dat$eaf_corrected = dat$observed_EAF - dat$shift

ggplot_list <- list()

for (name in sample_codes) {
  data.viz <- dat %>% filter(experimental_id == name)
  p1 <- ggplot(data.viz, aes(x = eaf_corrected)) +
    ggtitle(name) +
    geom_histogram() +
    geom_vline(xintercept = data.viz$eaf_max, color = "red") +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
          axis.title.x = element_text(size = 10),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(0.05, "inch"),
          legend.background = element_rect(fill = "transparent"))
  
  # Store the ggplot in the list
  ggplot_list[[name]] <- p1
}

fig3 = ggplot_list[[1]] + ggplot_list[[2]] + ggplot_list[[3]] + ggplot_list[[4]] + ggplot_list[[5]] + ggplot_list[[6]] 

fig3


#correct any remaining negative values so that they are equal to zero
dat[dat$eaf_corrected < 0, "eaf_corrected"] <- 0

# Loop through each value in sample_codes
results_df <- data.frame()

for (name in sample_codes) {
  dat3 <- dat %>% filter(experimental_id == name)
  replacement <- median(dat3$eaf_max - 0.02)
  dat4 <- dat3 %>% mutate(eaf_corrected = ifelse(eaf_corrected > eaf_max, replacement, eaf_corrected))
  results_df <- rbind(results_df, dat4)
}

#visualize
ggplot_list <- list()

for (name in sample_codes) {
  data.1 <- results_df %>% filter(experimental_id == name)
  p1 <- ggplot(data.1, aes(x = eaf_corrected)) +
    ggtitle(name) +
    geom_histogram() +
    geom_vline(xintercept = data.1$eaf_max, color = "red") +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
          axis.title.x = element_text(size = 10),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(0.05, "inch"),
          legend.background = element_rect(fill = "transparent"))
  ggplot_list[[name]] <- p1
}

fig4 = ggplot_list[[1]] + ggplot_list[[2]] + ggplot_list[[3]] + ggplot_list[[4]] + ggplot_list[[5]] + ggplot_list[[6]] 

fig4

write.csv(results_df, "16S/results/qsip_results_18O_outliers_removed.csv")
