## Tables for PC Relationships

library(tidyverse)
library(ggplot2)
library(phyloseq)
library(lmerTest)

rm(list=ls())
setwd("/data/Users/kmccauley/AZ_PPG/R21/Manuscript1/")
phy <- readRDS("/data/Users/kmccauley/AZ_PPG/R21/Manuscript1/analysis_phyloseq.rds")
grps <- read.csv("/data/Users/kmccauley/AZ_PPG/R21/Manuscript1/hh_newTRQ-IL6_groups_dividedatzero_forKatie_9-24-21.csv", na.strings="")

sample_data(phy) <- data.frame(sample_data(phy)) %>% 
  rownames_to_column("SampleID") %>% 
  left_join(grps, by="househld_ad") %>% 
  column_to_rownames("SampleID")

## Determine househld_ads with duplicates
#dups <- allphy@sam_data$househld_ad[duplicated(allphy@sam_data$househld_ad)]
#phy <- subset_samples(allphy, househld_ad %in% dups)
sample_data(phy)$trq_split <- sample_data(phy)$RatioQ1vQ234.MEturquoise > 0
sample_data(phy)$il6_split <- sample_data(phy)$z_IL6LPCC.y <= 0

table(sample_data(phy)$il6_split, sample_data(phy)$lohi6,useNA="a")
table(sample_data(phy)$newT6, sample_data(phy)$il6trq,useNA="a")

##Loop through variable and get p-values for each of the first three principal coordinates
myvars <- names(sample_data(phy))[!names(sample_data(phy)) %in% c("househld_ad", "AZ_ID", "UCSF_ID", "Axis.1","Axis.2","Axis.3","Axis.4","Axis.5")]
## New, split variables:

res <- list()
res2 <- list()
vartype <- list()
for(i in myvars) {
  if(length(unique(sample_data(phy)[,i])@.Data[[1]])>1) {
    for(j in 1:5) {
      print(i)
      myformula <- paste0("Axis.", j, " ~ ", i, " + age + (1|househld_ad)")
      res[[i]][[j]] <- anova(lmer(myformula, data=data.frame(sample_data(phy))))$`Pr(>F)`[[1]]
      myformula2 <- paste0("Axis.", j, " ~ ", i, "* age + (1|househld_ad)")
      res2[[i]][[j]] <- anova(lmer(myformula2, data=data.frame(sample_data(phy))))$`Pr(>F)`[3]
      vartype[[i]] <- class(data.frame(sample_data(phy))[,i])
    }
  }
}
wunifrac_results <- as.data.frame(do.call(rbind, res)) %>%
  rownames_to_column("Name of Variable") %>%
  mutate_at(vars(starts_with("V")), round, digits=3) %>% 
  rowwise() %>%
  filter(`Name of Variable` %in% c("collapsed_ratios","asthma2to9yr","domGenus","newT6","bigt","lohi6")) %>%
  rename("PC1"="V1", "PC2"="V2","PC3"="V3","PC4"="V4","PC5"="V5") %>%
  write.csv("wunifrac_results_PC1thru3_new.csv", row.names = FALSE)

wunifrac_results_int <- as.data.frame(do.call(rbind, res2)) %>%
  rownames_to_column("Name of Variable") %>%
  mutate_at(vars(starts_with("V")), round, digits=3) %>% 
  rowwise() %>%
  filter(`Name of Variable` %in% c("collapsed_ratios","asthma2to9yr","domGenus","newT6","bigt","lohi6")) %>%
  rename("PC1"="V1", "PC2"="V2","PC3"="V3","PC4"="V4","PC5"="V5") %>%
  write.csv("wunifrac_int_PC1thru3_new.csv", row.names = FALSE)
