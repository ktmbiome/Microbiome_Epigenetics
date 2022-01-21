## Analysis Phyloseq
library(tidyverse)
library(ggplot2)
library(lmerTest)
library(phyloseq)

rm(list=ls())

setwd("/data/Users/kmccauley/AZ_PPG/R21/Manuscript1/")
source("/data/Users/kmccauley/LabCode/DESeqFunctions.R")
phy <- readRDS("/data/Users/kmccauley/AZ_PPG/R21/multiply_rarefy/rarefied_otu_phyloseq.rds")
map_data <- read.csv("/data/Users/kmccauley/AZ_PPG/R21/SUBSET_methylMEs_IIS_microbiome_manifest_de-identified_2020-09-20.csv")
map_data2 <- map_data %>% 
  mutate_all(list(~na_if(., ""))) %>%
  mutate(cseason4 = case_when(
    collectionmonth %in%  c("September","October","November") ~ "Fall",
    collectionmonth %in%  c("December","January","February")  ~ "Winter",
    collectionmonth %in%  c("March","April","May")  ~ "Spring",
    TRUE ~ "Summer")) %>% 
  mutate(collapsed_ratios = recode(pmratioQ4_MASTN_ad2020, `0`="Q1",`1`="Q234",`2`="Q234",`3`="Q234")) %>% 
  mutate(cat_ratios = recode(pmratioQ4_MASTN_ad2020, `0`="Q1",`1`="Q2",`2`="Q3",`3`="Q4")) %>% 
  mutate(age = as.numeric(sub("mo","",VisitType))) %>%
  mutate(z_IL6LPCC = as.numeric(as.character(z_IL6LPCC))) %>% 
  mutate(z_TNFLPCC = as.numeric(as.character(z_TNFLPCC))) %>% 
  mutate(z_li1blpcc = as.numeric(as.character(z_li1blpcc))) %>% 
  mutate(lifncrat3m = as.numeric(as.character(lifncrat3m))) %>%
  mutate(lifncrat3m_q = ntile(lifncrat3m, 4)) %>%
  mutate_all(~na_if(., "")) %>% 
  mutate_all(~na_if(., "#N/A")) %>% 
  group_by(househld_ad) %>% 
  arrange(age) %>% 
  mutate(sampcount=row_number()) %>% 
  as.data.frame() %>% 
  column_to_rownames("UCSF_ID")

## Made the astute observation that we definitely have repeated measures here.
## CpGs were calculated  AT BIRTH (so yes, there's only one timepoint for that data). Therefore, we are asking a development question with the nasal data.

phy_az <- phy
sample_data(phy_az) <- sample_data(map_data2)

sample_data(phy_az)$il6trq_ex <- sample_data(phy_az)$il6trq
sample_data(phy_az)$il6trq_ex[!sample_data(phy_az)$il6trq %in% c("hi6-lowT", "low6-hiT")] <- NA

phy.nasal <- subset_samples(phy_az, Material %in% "NarSwb")
phy.nasal.filt <- subset_taxa(phy.nasal, taxa_sums(phy.nasal)>0)
phy_tree(phy.nasal.filt) <- ape::root(phy_tree(phy.nasal.filt), 1, resolve.root=TRUE)
phy.nasal.filt

final <- data.frame(tax_table(phy.nasal.filt)@.Data) %>%
  rownames_to_column("SVno") %>% 
  mutate_all(as.character) %>%
  column_to_rownames("SVno")
final$Genus[final$Genus %in% "NA"] <- final$Family[final$Genus %in% "NA"]
final$Genus[final$Genus %in% "NA"] <- final$Order[final$Genus %in% "NA"]
final$Genus[final$Genus %in% "NA"] <- final$Class[final$Genus %in% "NA"]
tax_table(phy.nasal.filt) <- as.matrix(final)

phy.genus <- tax_glom(phy.nasal.filt, "Genus")
sample_data(phy.nasal.filt)$domGenus <- c(tax_table(phy.genus)[apply(otu_table(phy.genus),2, which.max)][,"Genus"])
other_tax <- names(table(sample_data(phy.nasal.filt)$domGenus)[table(sample_data(phy.nasal.filt)$domGenus) < 2])
sample_data(phy.nasal.filt)$domGenus_other <- sample_data(phy.nasal.filt)$domGenus
sample_data(phy.nasal.filt)$domGenus_other[sample_data(phy.nasal.filt)$domGenus_other %in% other_tax] <- "Other"
sample_data(phy.nasal.filt)$domGenus_pct <- apply(otu_table(phy.genus),2, max)/apply(otu_table(phy.genus),2, sum)

wunidist <- as.matrix(phyloseq::distance(phy.nasal.filt, method="wunifrac"))
dist.for.lme <- ape::pcoa(wunidist)$vectors[,c(1:5)]
int.dat <- merge(sample_data(phy.nasal.filt), dist.for.lme, by=0)
row.names(int.dat) <- int.dat$Row.names
int.dat$Row.names <- NULL
sample_data(phy.nasal.filt) <- sample_data(int.dat)

#write.csv(int.dat, file="SUBSET_with_MicrobiomeVariables_2020-09-11.csv")

phy.nasal.filt
saveRDS(phy.nasal.filt, "analysis_phyloseq.rds")
