## Figure Related to Relationship Between Microbiome and Methylation Patterns


library(tidyverse)
library(ggplot2)
library(lmerTest)

rm(list=ls())


setwd("/data/Users/kmccauley/AZ_PPG/R21/Manuscript1/")
source("/data/Users/kmccauley/LabCode/DESeqFunctions.R")
phy <- readRDS("/data/Users/kmccauley/AZ_PPG/R21/Manuscript1/analysis_phyloseq.rds")

merged.data <- data.frame(sample_data(phy))
merged.data <- merged.data[merged.data$age >= 24, ]
for(i in c(24,36)) {
  ss.data <- merged.data[merged.data$age == i, ]
  print(summary(lm(Axis.1 ~ RatioQ1vQ234.MEturquoise + collapsed_ratios, data=ss.data))$coef[2,4])
}

merged.data$VisitType <- factor(merged.data$VisitType, levels=c("24mo","36mo"))

lab_dat <- data.frame(label=c("P=0.020","P=0.583"), VisitType=c("24mo","36mo"))

partA <- ggplot(merged.data, aes(y=RatioQ1vQ234.MEturquoise, x=Axis.1)) +
  xlab("Weighted UniFrac PC1") +
  geom_text(data=lab_dat,aes(x=0, y=0.23, label=label)) +
  ylab("TRQ Eigengene") +
  geom_point() + 
  facet_wrap(~VisitType) +
  geom_smooth(method="lm", linetype=2, color="black")
partA


ggsave("ModuleEigengenePlot_new.pdf", partA, device=cairo_pdf, height=4, width=6)
