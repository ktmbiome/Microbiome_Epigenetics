## Figure 1. Nasal microbiome development over the first three years of life.

library(tidyverse)
library(ggplot2)
library(lmerTest)
library(phyloseq)

rm(list=ls())
theme_set(theme_classic())
ratio_colors = c("grey","black")
trq_colors = c("grey","black")

setwd("/data/Users/kmccauley/AZ_PPG/R21/Manuscript1/")
source("/data/Users/kmccauley/LabCode/SigTaxa/DESeqFunctions.R")
phy <- readRDS("/data/Users/kmccauley/AZ_PPG/R21/Manuscript1/analysis_phyloseq.rds")
grps <- read.csv("/data/Users/kmccauley/AZ_PPG/R21/Manuscript1/hh_newTRQ-IL6_groups_dividedatzero_forKatie_9-24-21.csv", na.strings="")

sample_data(phy) <- data.frame(sample_data(phy)) %>% 
  rownames_to_column("SampleID") %>% 
  left_join(grps, by="househld_ad") %>% 
  column_to_rownames("SampleID")

sample_data(phy)$il6trq_ex <- factor(sample_data(phy)$il6trq_ex, labels=c("High IL-6 / Low TRQ", "Low IL-6 / High TRQ"))
sample_data(phy)$trq_split <- sample_data(phy)$RatioQ1vQ234.MEturquoise >= 0
sample_data(phy)$il6_split <- sample_data(phy)$z_IL6LPCC.x >= 0
sample_data(phy)$il6trq_onevall <- sample_data(phy)$newT6 %in% "hi6,lowT"
sample_data(phy)$il6trq_extremes <- sample_data(phy)$newT6
sample_data(phy)$il6trq_extremes[!sample_data(phy)$il6trq_extremes %in% c("hi6,lowT","low6,hiT")] <- NA

make_interactionPlots <- function(data, variable, axis.num, colors, note=variable, mylabs=waiver(), lines=c("solid","solid")) {
  myfunction <- paste0("Axis.", axis.num, " ~ ", variable, "* age + (1|househld_ad)")
  mydata <- data[!is.na(data[,variable]),]
  print(dim(mydata))
  plot.model <- lmer(as.formula(myfunction), data=mydata)
  aov <- anova(plot.model)
  print(aov)
  mydata$predicted <- predict(plot.model)

actual.data <- ggplot(mydata, aes_string(x="age", y=paste0("Axis.", axis.num), fill=variable, color=variable, linetype=variable)) +
  ylab(paste0("Weighted UniFrac PC", axis.num)) +
  annotate("text", x=Inf, y=Inf, hjust=1.3, vjust=1.3, label=paste0("Pint = ", round(aov$`Pr(>F)`[[nrow(aov)]],3))) +
  geom_smooth(aes(y=predicted), method="lm", se=FALSE) + 
  geom_point(shape=21, size=3, color="black") +
  xlab("Child's Age at Collection") +
  scale_fill_manual(note, values=colors, labels=mylabs) + 
  scale_color_manual(note, values=colors, labels=mylabs) +
  scale_linetype_manual(note, values=lines, labels=mylabs)
return(actual.data)
}

anly.dat <- data.frame(sample_data(phy))
coords <- data.frame(xmin=1, xmax=Inf, ymin=c(-0.2), ymax=c(0))
myfig1 <- make_interactionPlots(anly.dat, "collapsed_ratios", 3, colors=ratio_colors, "Maternal Ratio")
fig1_colorized <- myfig1 +
  annotate("rect", xmin=1, xmax=Inf, ymin=-0.2, ymax=0, alpha=0.5, fill="darkblue") +
  annotate("rect", xmin=1, xmax=Inf, ymin=0, ymax=0.5, alpha=0.5, fill="goldenrod")
## Reorganize the layers so that the rectangles are in the background
fig1_colorized$layers <- c(fig1_colorized$layers[[4]], fig1_colorized$layers[[5]], fig1_colorized$layers[[1]], fig1_colorized$layers[[2]], fig1_colorized$layers[[3]])
fig1_colorized



## the goal of these figures is to show significant interactions, not main effects
anly.dat$wza29 <- factor(anly.dat$wza29, labels=c("no","yes"))
myfig_wza <- make_interactionPlots(anly.dat, "wza29", 2)
myfig_il6ex <- make_interactionPlots(anly.dat, "il6trq_extremes", 2, colors=trq_colors, "IL6/Turquoise Ratio")
anly.dat$z_IL6LPCC_q <- factor(ntile(anly.dat$z_IL6LPCC, 2), labels=paste0("Q", 1:2))
myfig_il6c <- make_interactionPlots(anly.dat, "z_IL6LPCC_q", 1, trq_colors)
anly.dat$RatioQ1vQ234.MEturquoise_q <- factor(ntile(anly.dat$RatioQ1vQ234.MEturquoise, 4), labels=paste0("Q", 1:4))
myfig_ratq <- make_interactionPlots(anly.dat, "RatioQ1vQ234.MEturquoise_q", 3, c("red","green","blue","purple"))
make_interactionPlots(anly.dat, "cseason4", 1)
myfig_asth <- make_interactionPlots(anly.dat, "asthma2to9yr", 2, trq_colors, "Childhood Asthma", c("No","Yes"))
anly.dat$bigt <- factor(anly.dat$bigt)
anly.dat$lohi6 <- factor(anly.dat$lohi6)
myfig_trq0 <- make_interactionPlots(anly.dat, "bigt", 3, trq_colors, "Turquoise Module", mylabs=c("Negative Eigenvector","Positive Eigenvector"))
myfig_il60 <- make_interactionPlots(anly.dat, "lohi6", 2, trq_colors, "IL-6 at Birth", mylabs=c("Positive Z-score","Negative Z-score"))

myfig_t6 <- make_interactionPlots(anly.dat, "newT6", 2, c(trq_colors, trq_colors), "IL6/TRQ Groups", lines=c("solid","solid","dashed","dashed"))
myfig_t6

myfig_t6_v2 <- make_interactionPlots(anly.dat, "il6trq_onevall", 2,trq_colors, "IL6/TRQ Groups", mylabs=c("All Others","hi6,lowT"))

# ggsave("IntEffect_Ratio.pdf", myfig1, device="pdf", height=6, width=6)
# ggsave("IntEffect_TRQmed.pdf", myfig_trq0, device="pdf", height=6, width=6)
# ggsave("IntEffect_IL6med.pdf", myfig_il60, device="pdf", height=6, width=6)
# ggsave("IntEffect_Asthma.pdf", myfig_asth, device="pdf", height=6, width=6)
# ggsave("IntEffect_Newil6trq.pdf", myfig_t6, device="pdf", height=6, width=6)
# ggsave("IntEffect_newT6_vAll.pdf", myfig_t6_v2, device="pdf",height=6,width=6)
# ggsave("IntEffect_il6trq_EXTR.pdf",myfig_il6ex, device="pdf", height=6, width=6)

make_effectPlots <- function(data, variable, axis.num, colors, note=variable, mylabs=waiver()) {
  myfunction <- paste0("Axis.", axis.num, " ~ ", variable, " + age + (1|househld_ad)")
  mydata <- data[!is.na(data[,variable]),]
  plot.model <- lmer(as.formula(myfunction), data=mydata)
  print(summary(plot.model)$coef)
  mydata$predicted <- predict(plot.model)
  
  actual.data <- ggplot(mydata, aes_string(x="age", y=paste0("Axis.", axis.num), fill=variable, color=variable)) +
    ylab(paste0("Weighted UniFrac PC", axis.num)) +
    annotate("text", x=Inf, y=Inf, hjust=1, vjust=1.3, label=paste0("P = ", round(summary(plot.model)$coef[2,5],3))) +
    geom_smooth(aes(y=predicted), method="lm", se=FALSE) + 
    geom_point(shape=21, size=3, color="black") +
    xlab("Child's Age at Sample Collection") +
    scale_fill_manual(note, values=colors, labels=mylabs) + 
    scale_color_manual(note, values=colors, labels=mylabs)
  return(actual.data)
}
myfig_main_rat <- make_effectPlots(anly.dat, "collapsed_ratios", 2, trq_colors, note="Maternal Ratio")
myfig_main_trq <- make_effectPlots(anly.dat, "trq_split", 2, c(0="grey",0="black"), note="TRQ (median)", mylabs=c("Below Median","Above Median"))
myfig_main_il6 <- make_effectPlots(anly.dat, "il6_split", 2, trq_colors)
myfig_main_asth <- make_effectPlots(anly.dat, "asthma2to9yr", 2, trq_colors, note="Childhood Asthma")
myfig_main_newT6 <- make_effectPlots(anly.dat, "newT6", 2, c("#333333","#606060","#c0c0c0","#f2f2f2"), note="IL6/TRQ Groups")


myfig_main_trq <- make_effectPlots(anly.dat, "bigt", 2, trq_colors, "Turquoise Module", mylabs=c("Negative Eigenvector","Positive Eigenvector"))
myfig_main_il6 <- make_effectPlots(anly.dat, "lohi6", 2, trq_colors, "IL-6 at Birth", mylabs=c("Positive Z-score","Negative Z-score"))

#ggsave("MainEffect_Ratio.pdf", myfig_main_rat, device="pdf", height=6, width=6)
#ggsave("MainEffect_TRQ_zero.pdf", myfig_main_trq, device="pdf", height=6, width=6)
#ggsave("MainEffect_IL6_zero.pdf", myfig_main_il6, device="pdf", height=6, width=6)
#ggsave("MainEffect_Asthma.pdf", myfig_main_asth, device="pdf", height=6, width=6)
#ggsave("MainEffect_Newil6trq.pdf", myfig_main_newT6, device="pdf", height=6, width=6)

library(ggpubr)
row1 <- ggarrange(NULL, myfig_main_trq, myfig_main_asth, NULL,nrow=1, widths=c(1,2,2,1), labels=c("","A","B",""))
row2 <- ggarrange(myfig1, myfig_trq0, myfig_il60, nrow=1, labels=c("C","D","E"))
final_plot <- ggarrange(row1, row2, nrow=2)
ggsave("Main_Int_Effect_Figure.pdf", final_plot, width = 15, height=7, device=cairo_pdf)



#### Make PC Plots
merged.otus <- merge(t(otu_table(phy)), sample_data(phy), by=0)
otus <- rownames(otu_table(phy))
sav <- list()
for(i in 1:length(otus)) {
  sav1 <- otus[i]
  sav2 <- cor(merged.otus[,paste0("Axis.1")], merged.otus[,otus[i]], method="spearman")
  sav3 <- cor(merged.otus[,paste0("Axis.2")], merged.otus[,otus[i]], method="spearman")
  sav4 <- cor(merged.otus[,paste0("Axis.3")], merged.otus[,otus[i]], method="spearman")
  sav[[i]] <- cbind(sav1, sav2, sav3, sav4)
}
library(tidyverse)

cor.results <- data.frame(do.call(rbind, sav)) %>%
  mutate(sav2 = as.numeric(as.character(sav2))) %>% 
  mutate(sav3 = as.numeric(as.character(sav3))) %>% 
  mutate(sav4 = as.numeric(as.character(sav4)))
final <- merge(cor.results, tax_table(phy)@.Data, by.x="sav1", by.y=0) %>% 
  mutate_at(c("Genus","Family","Order","Class"), as.character) %>% 
  data.frame()
final$Genus[final$Genus %in% NA] <- final$Family[final$Genus %in% NA]
final$Genus[final$Genus %in% NA] <- final$Order[final$Genus %in% NA]
final$Genus[final$Genus %in% NA] <- final$Class[final$Genus %in% NA]
final$Genus[final$Genus %in% NA] <- final$Phylum[final$Genus %in% NA]
tail(final, 30)

final$labels <- paste0(final$Genus, " (", final$sav1, ")")

mylabs <- function(mytext) {
  sapply(
  strsplit(as.character(mytext), " "), 
  function(x) parse(text = paste0("italic('", x[1], "')~", x[2]))
)}

##labs(final$labels[reorder(final$labels[abs(final$sav2) > 0.25], sav2)])

# ggplot(data, aes(name, value)) + 
#   geom_col() + coord_flip() +
#   scale_x_discrete(labels = labs)

PC1 <- arrange(subset(final, abs(sav2) > 0.25), sav2)
PC1$labels <- factor(PC1$labels, levels=PC1$labels)
figPC1 <- ggplot(PC1, aes(x=sav2, y=labels, fill=sav2 > 0)) + 
  geom_bar(stat="identity") +
  geom_vline(xintercept=c(0.5, -0.5), linetype="dashed") +
  scale_fill_manual(values=c("darkblue","goldenrod")) +
  ylab("") +
  xlab("Spearman Correlation with PC1") +
  guides(fill=FALSE) +
  expand_limits(x=c(-1,1)) +
  scale_y_discrete(labels=mylabs(levels(PC1$labels)))
figPC1

figPC2 <- ggplot(subset(final, abs(sav3) > 0.25), aes(x=sav3, y=reorder(paste0(Genus, " (", sav1, ")"), sav3), fill=sav3 > 0)) + 
  geom_bar(stat="identity") +
  geom_vline(xintercept=c(0.5, -0.5), linetype="dashed") +
  scale_fill_manual(values=c("goldenrod")) +
  ylab("") +
  xlab("Spearman Correlation with PC2") +
  guides(fill=FALSE) +
  expand_limits(x=c(-1,1))

PC3 <- arrange(subset(final, abs(sav4) > 0.25), sav4)
PC3$labels <- factor(PC3$labels, levels=PC3$labels)
figPC3 <- ggplot(PC3, aes(x=sav4, y=labels, fill=sav4 > 0)) + 
  geom_bar(stat="identity") +
  geom_vline(xintercept=c(0.5, -0.5), linetype="dashed") +
  scale_fill_manual(values=c("darkblue","goldenrod")) +
  ylab("") +
  xlab("Spearman Correlation with PC3") +
  guides(fill=FALSE) +
  expand_limits(x=c(-1,1)) +
  scale_y_discrete(labels=mylabs(levels(PC3$labels)))

# ggsave("PCLoadings_PC123.eps", gridExtra::grid.arrange(figPC1, figPC2, figPC3, nrow=1), device="eps", height=5, width=15, dpi = 300)
# ggsave("PCLoadings_PC1.pdf", figPC1, device=cairo_pdf, height=5, width=6)
# ggsave("PCLoadings_PC3.pdf", figPC3, device=cairo_pdf, height=5, width=6)
##ggsave("Figure_1.pdf", myfig1, device=cairo_pdf, height=5, width=6)

## Make Panel C
merged.data <- data.frame(sample_data(phy))
merged.data <- merged.data[merged.data$age >= 24, ]
for(i in c(24,36)) {
  ss.data <- merged.data[merged.data$age == i, ]
  print(summary(lm(Axis.1 ~ RatioQ1vQ234.MEturquoise + collapsed_ratios, data=ss.data))$coef[2,4])
}

merged.data$VisitType <- factor(merged.data$VisitType, levels=c("24mo","36mo"))

lab_dat <- data.frame(label=c("P=0.020","P=0.583"), VisitType=c("24mo","36mo"))

trq_PC1_corr <- ggplot(merged.data, aes(y=RatioQ1vQ234.MEturquoise, x=Axis.1)) +
  xlab("Weighted UniFrac PC1") +
  geom_text(data=lab_dat,aes(x=0, y=0.23, label=label)) +
  ylab("TRQ Eigengene") +
  geom_point() + 
  facet_wrap(~VisitType) +
  geom_smooth(method="lm", linetype=2, color="black")
trq_PC1_corr


library(ggpubr)
row1 <- ggarrange(fig1_colorized, figPC3, labels=c("A","B"))
row2 <- ggarrange(trq_PC1_corr, figPC1, labels=c("C","D"), widths=c(1,0.7))
final_plot <- ggarrange(row1, row2, nrow=2)
ggsave("Figure4_combined.pdf", final_plot, device=cairo_pdf, height=7, width=10)
