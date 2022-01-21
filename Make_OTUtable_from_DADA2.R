# DADA2 Post-Processing.
#Takes the files created in DADA2 and makes them useful/qiime compatible. The names of these files are consistent with the names given in the output files from my "DADA2_Script.R" code.
# I also rename sequences as "OTUs" here, and save the sequences elsewhere. 
# I suggest running this with "nohup", since the "optimal tree" generation can take up to 12 hours. However, the OTU table and a temporary tree should print out relatively quickly.

# Author: Katie McCauley

rm(list=ls())
setwd(".")
dada2_otu <- read.table("otutable.txt", header=TRUE, check.names=F, comment="", sep="\t", row.names=1)
dada2_tax <- read.table("tax_table.txt", header=TRUE, check.names=F, comment="", sep="\t") ## Moving forward with DECIPHER because the DADA2 taxonomy didn't generate for about 2000 sequences.
dada2_seqs <- read.table("seqs.txt", header=F, comment="", sep="\t")

#Make the sample names actual sample names instead of file names:
listnames <- strsplit(names(dada2_otu), split="_", fixed=TRUE)
names(dada2_otu) <- sapply(lapply(listnames, function(x) x[c(1) ]), paste, collapse="_") ## Maintains the run information in the sample name so that we can figure out duplicates

## Identify duplicate samples
dupnames <- data.frame(do.call(rbind, strsplit(names(dada2_otu), split="_", fixed=TRUE)))
dupnames$mydups <- dupnames[,1] %in% dupnames[,1][duplicated(dupnames[,1])]
rownames(dupnames) <- names(dada2_otu)
## How can I get sample-specific data from 

## Move to Phyloseq structure...
library(phyloseq)
phy <- phyloseq(otu_table(dada2_otu, taxa_are_rows=TRUE), sample_data(dupnames), tax_table(as.matrix(dada2_tax)))
sequences <- Biostrings::DNAStringSet(taxa_names(phy))
names(sequences) <- taxa_names(phy)
phy <- merge_phyloseq(phy, sequences)
phy
taxa_names(phy) <- paste0("SV_", 1:length(taxa_sums(phy)))

phy_filtB <- subset_taxa(phy, Kingdom %in% "Bacteria")
phy_filt <- prune_taxa(taxa_sums(phy_filtB) > 0.000001*sum(taxa_sums(phy_filtB)), phy_filtB)
phy_filt
saveRDS(phy_filt, file="az_dada2_phy_init.rds")

library(phangorn)
library(msa)
library(DECIPHER)
phy_filt_tree2 <- phy_filt_tree <- phy_filt
print("Running Alignment")
alignment <- AlignSeqs(DNAStringSet(refseq(phy_filt_tree)), anchor=NA)
print("Running phangorn")
phang.align <- as.phyDat(as(alignment, "matrix"), type="DNA")
print("Making distance matrix")
dm <- dist.ml(phang.align)
print("Making tree from distance matrix")
treeNJ <- NJ(dm)
print("Fitting")
fit <- pml(treeNJ, data=phang.align)
phy_filt_tree <- merge_phyloseq(phy_filt, phy_tree(fit$tree))
saveRDS(phy_filt_tree, file="az_dada2_phy_tree1.rds")
fitGTR <- update(fit, k=4, inv=0)
print("Optimizing Fit")
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
          rearrangement = "stochastic", control = pml.control(trace = 0))
write.tree(fitGTR$tree, "dada2_optim_tree.tre")
phy_filt_tree <- merge_phyloseq(phy_filt, phy_tree(fitGTR$tree))
saveRDS(phy_filt_tree, file="az_dada2_phy_final.rds")
