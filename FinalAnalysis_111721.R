## Final analysis with specific CpGs and specific taxa.

library(tidyverse)
library(ggplot2)
library(lmerTest)
library(phyloseq)

rm(list=ls())
theme_set(theme_classic())

setwd("/data/Users/kmccauley/AZ_PPG/R21/Manuscript1/")
phy <- readRDS("/data/Users/kmccauley/AZ_PPG/R21/Manuscript1/analysis_phyloseq.rds")
dmcs <- read.csv("/data/Users/kmccauley/AZ_PPG/R21/asthmaDMCs_TRQ_matrix_forKatie.csv", row.names=1)
keeps <- read.csv("/data/Users/kmccauley/AZ_PPG/R21/trq_asthmaDMCs_promoter-enhancers_Nov9.csv", row.names=1)
taxa.of.int <- c("SV_1","SV_196","SV_197","SV_247","SV_329","SV_6","SV_2","SV_23","SV_4")
phy.filt <- prune_taxa(taxa = taxa.of.int, phy)
phy24 <- subset_samples(phy.filt, VisitType %in% "24mo")
table(sample_data(phy.filt)$age)
phy.data <- merge(t(otu_table(phy24)), sample_data(phy24), by=0)
comb <- merge(phy.data, dmcs, by.x="househld_ad", by.y=0)


keeps$Gene_m <- as.character(keeps$geneSymbol)
keeps$Gene_m[keeps$Gene_m %in% ""] <- as.character(keeps$ProbeID)[keeps$Gene_m %in% ""]
dim(comb)
cpg.names <- names(dmcs)

res <- list()
pval <- list()
for(i in taxa.of.int) {
  for(j in cpg.names) {
    res[[i]][[j]] <- cor.test(comb[,i], comb[,j], method="s", exact=FALSE)$estimate
    pval[[i]][[j]] <- cor.test(comb[,i], comb[,j], method="s", exact=FALSE)$p.value
  }
}
pval.df <- data.frame(pval)
res.df <- data.frame(res)

res.df.filt <- res.df[rownames(res.df) %in% keeps$Gene_m & rowSums(pval.df < 0.05) > 0, colSums(pval.df < 0.05) > 0]
pval.df.filt <- round(pval.df[rownames(res.df) %in% keeps$Gene_m & rowSums(pval.df < 0.05) > 0, colSums(pval.df < 0.05) > 0], 2)
pval.df.filt[pval.df.filt > 0.05] <- ""
## Instead of showing p-values, maybe show an asterisk for significant correlations
#pval.df.filt <- pval.df[rowSums(pval.df < 0.05) > 0, colSums(pval.df < 0.05) > 0]
#pval.df.filt1 <- pval.df.filt
#pval.df.filt1[pval.df.filt1 < 0.05] <- "o"
#pval.df.filt1[!pval.df.filt1 == "o"] <- ""

colannos <- data.frame(SVs=names(res.df.filt), PCDirection=c(rep("Negative",4), rep("Positive",4)))
rownames(colannos) <- colannos$SVs
colannos$SVs <- NULL

colors=list(
  PCDirection=c("Positive"="goldenrod", "Negative"="darkblue")
)


pheatmap::pheatmap(res.df.filt,
                   clustering_method="complete",
                   display_numbers=pval.df.filt, 
                   fontsize_number=7,
                   annotation_col = colannos,
                   #annotation_row = row_annos,
                   annotation_colors=colors,
                   border_color="white",
                   filename = "CpG_SV_Correlation.pdf",
                   height = 8, width=7)
library(ComplexHeatmap)
library(RColorBrewer)
ha=HeatmapAnnotation(df=data.frame(colannos),col=list(PCDirection=c("Positive"="goldenrod","Negative"="darkblue")))
mylabels <- c("Moraxella (SV_1)","Moraxella (SV_196)","Moraxella (SV_197)","Moraxella (SV_247)","Corynebacterium (SV_6)","Dolosigranulum (SV_2)", "Acinetobacter (SV_23)","Staphylococcus (SV_4)")
pdf("CpG_SV_Correlation_plot_111721.pdf", height=10, width=7)
Heatmap(as.matrix(res.df.filt), col=rev(brewer.pal(n = 7, name = "RdYlBu")),
        clustering_method_rows="complete", 
        clustering_method_columns = "complete", 
        name="Correlation\nCoefficient",
        cell_fun=function(j,i,x,y,width,height,fill) {
          grid.text(sprintf("%s", pval.df.filt[i,j]),x,y,gp=gpar(fontsize=10))
        },
        top_annotation=ha,
        rect_gp=gpar(col="white", lwd=1),
        column_names_max_height = max_text_width(mylabels),
        column_labels=mylabels)
dev.off()
write.csv(pval.df.filt, "P_value_SV_by_CpG.csv")
write.csv(pval.df.filt1, "Sig_SV_by_CpG.csv")

## What proportion of my CpGs are MRE's?
row_anno_dat.ss <- row_anno_dat[rownames(res.df.filt),]
sum(!is.na(row_anno_dat.ss$MR))/nrow(row_anno_dat.ss)
