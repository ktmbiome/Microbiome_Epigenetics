## Make Alpha Rarefaction Curves

library(vegan)
library(phyloseq)

set.seed(123)
phy <- readRDS("/data/Users/kmccauley/AZ_PPG/dada2_output/phyloseq_noneg.rds")
phy <- subset_samples(phy, sample_sums(phy) < 10000)
arare <- vegan::rarecurve(as.matrix(t(otu_table(phy))), step=500, cex=0.5)

samplenames <- unlist(lapply(arare, function(x) last(names(attr(x, "Subsample")))))
print(arare)
plot(arare)
### 3000 looks good everythhing that's lower abundance evens out around there. No additional information is gained.
