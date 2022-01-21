# Code to remove negative controls.
# Katie McCauley
# Is it possible to develop this code so that it plays well with the pipeline so far.
#Maybe use the input CSV to determine which samples belong to each run and do run-specific negative control filtering?

library(phyloseq)

setwd("/data/Users/kmccauley/AZ_PPG/dada2_output/")
phy <- readRDS("az_dada2_phy_final.rds")

neg.dat <- subset_samples(phy, grepl("NTC|EMPTY", sample_names(phy)))
samp.dat <- subset_samples(phy,!grepl("NTC|EMPTY", sample_names(phy)))
negsA <- rowSums(otu_table(neg.dat)>0)/ncol(otu_table(neg.dat))

#The proportion of samples with counts for each OTU
sampsA <- rowSums(otu_table(samp.dat)>0)/ncol(otu_table(samp.dat))
dat <- merge(negsA, sampsA, by=0)
names(dat) <- c("#OTU ID","negs","samps")
dat2 <- merge(dat, tax_table(phy)@.Data, by.x="#OTU ID",by.y=0)
levels(dat2$Genus)[levels(dat2$Genus) %in% "NA"] <- "Other" #Most-frequent OTU classifications

library(ggplot2)
pre.cleaning <- ggplot(dat2, aes(negs,samps)) + 
  geom_point() + 
  ylab("Proportion of Samples") + 
  xlab("Proportion of Negative Controls") + 
  ggtitle("") + 
  ggrepel::geom_label_repel(data=subset(dat2, negs > 0.2 & samps < 0.6), aes(label=Genus)) +
  labs(tag="A")
#ggsave("Pre_Cleaning_Figure.pdf", pre.cleaning, device="pdf", height=6, width=7)

outright.drop <- dat2$`#OTU ID`[dat2$negs > 0.20 & dat2$samps < 0.60]

#ord = outright.drop
otu.ord <- subset_taxa(phy, !taxa_names(phy) %in% outright.drop)

#otu.sub
neg.dat2 <- subset_taxa(neg.dat, taxa_sums(neg.dat)>0 & !taxa_names(neg.dat) %in% outright.drop)
mean.in.NTC <- ceiling(apply(otu_table(neg.dat2),1, mean)) # Want a round number, so I am taking the ceiling.
max.in.NTC <- apply(otu_table(neg.dat2),1, max)

otu.notax.mean <- t(otu_table(otu.ord))
for(i in names(mean.in.NTC)) {
  otu.notax.mean[, i] <- otu.notax.mean[, i] - mean.in.NTC[i]
  otu.notax.mean[, i][otu.notax.mean[, i] < 0] <- 0
}
otu.mean <- merge_phyloseq(otu.ord, otu_table(otu.notax.mean, taxa_are_rows=FALSE))

otu.notax.max <- t(otu_table(otu.ord))
for(i in names(max.in.NTC)) {
  otu.notax.max[, i] <- otu.notax.max[, i] - max.in.NTC[i]
  otu.notax.max[, i][otu.notax.max[,i] < 0] <- 0
}
otu.max <- merge_phyloseq(otu.ord, otu_table(otu.notax.max, taxa_are_rows=FALSE))


neg.dat <- subset_samples(otu.mean, grepl("NTC|EMPTY", sample_names(otu.mean)))
samp.dat <- subset_samples(otu.mean,!grepl("NTC|EMPTY", sample_names(otu.mean)))
negsA <- rowSums(otu_table(neg.dat)>0)/ncol(otu_table(neg.dat))
sampsA <- rowSums(otu_table(samp.dat)>0)/ncol(otu_table(samp.dat))
dat <- merge(negsA, sampsA, by=0)
names(dat) <- c("#OTU ID","negs","samps")
dat2 <- merge(dat, tax_table(phy)@.Data, by.x="#OTU ID",by.y=0)
levels(dat2$Genus)[levels(dat2$Genus) %in% "NA"] <- "Other" #Most-frequent OTU classifications

post.clean <- ggplot(dat2, aes(negs,samps)) + 
  geom_point() + 
  ylab("Proportion of Samples") + 
  xlab("Proportion of Negative Controls") + 
  ggtitle("") + 
  ggrepel::geom_label_repel(data=subset(dat2, negs > 0.05), aes(label=paste(Genus, `#OTU ID`))) +
  scale_color_brewer(palette="Paired") +
  labs(tag="B")
library(gridExtra)
ggsave("Combined_Cleaning_Figure.pdf", grid.arrange(pre.cleaning, post.clean, nrow=1), device="pdf", height=5, width=14)

saveRDS(samp.dat, "phyloseq_noneg.rds")
