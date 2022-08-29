# This script plots a heatmap with the ARG class prevalence per taxonomic phylum

# packages
library(data.table)
library(ggplot2)
library(dplyr)
library(viridis)
library(pheatmap)
library(stringr)

# load data
# load gOTU data, only want to analyze the MAGs that got clustered
gotu.final <- fread("data\\gotu_final.tsv", sep = "\t")
# only ANI 95
gotu.final <- gotu.final[which(gotu.final$ani == 95), ]

# load ARG data
deeparg <- fread("data\\deeparg.tsv", sep = "\t")
# only consider ARGs with 35% or higher percent identity
deeparg <- deeparg[deeparg$identity >= 35, ]

# load best bins to assign taxonomy to clusters
gotu.best <- fread("data\\gotu_best.tsv", sep = "\t")

# make sure taxonomy is consistent per cluster
taxonomy <- merge(gotu.final[, c("bin", "new_group")], gotu.best[, c("new_group", "taxonomy")], by = "new_group", all.x = T)

# extract phylum
taxonomy$phylum <- str_extract(taxonomy$taxonomy, "(?<=p__).*?(?=;)")
taxonomy$phylum <- gsub("_([[:alpha:]]{1,2})", "", taxonomy$phylum)

# put all data together in a table
arg.class.tax <- merge(gotu.final[, c("bin", "new_group")], taxonomy[, c("bin", "phylum")], by = "bin", all.x = T)
arg.class.tax <- merge(arg.class.tax, deeparg[, c("mag", "arg", "predicted.arg.class")], by.x = "bin", by.y = "mag", all.x = T)

# summarize data
arg.class.tax.cnt <- 
  merge(
    as.data.frame(arg.class.tax %>% 
      group_by(phylum, predicted.arg.class) %>% 
      summarize(
        gotu.arg.class.count = n_distinct(new_group, na.rm = T))
      ),
    as.data.frame(arg.class.tax %>% 
      group_by(phylum) %>% 
      summarize(
        gotu.count = n_distinct(new_group, na.rm = T))
      ),
    by = "phylum"
    
  )
arg.class.tax.cnt$prevalence <- arg.class.tax.cnt$gotu.arg.class.count/arg.class.tax.cnt$gotu.count
arg.class.tax.cnt <- arg.class.tax.cnt[!is.na(arg.class.tax.cnt$predicted.arg.class),]
# 
ac.tx <- matrix(0, 
                nrow = length(unique(arg.class.tax.cnt$predicted.arg.class)), 
                ncol = length(unique(arg.class.tax.cnt$phylum)))
rownames(ac.tx) <- unique(arg.class.tax.cnt$predicted.arg.class)
colnames(ac.tx) <- unique(arg.class.tax.cnt$phylum)

for (i in 1:nrow(ac.tx)) {
  for (j in 1:ncol(ac.tx)) {
    prev <- arg.class.tax.cnt$prevalence[
      which((arg.class.tax.cnt$predicted.arg.class == rownames(ac.tx)[i]) & (arg.class.tax.cnt$phylum == colnames(ac.tx)[j]))]
    if (length(prev) > 0) {
      ac.tx[i, j] <- prev
    }
  }
}

# color scale breaks
breaksList = seq(0, 1, by = 0.01)
# count number of mags in phylum
cnt.taxa <- as.data.frame(
  arg.class.tax %>% group_by(phylum) %>% summarize(unique_gotus = n_distinct(new_group))
)

# Modify ordering of the clusters using clustering callback option
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

# plot the heatmap
pdf(file = paste0("figures\\prevalence_argclass_phylum.pdf"),
    height = 11.69, width = 8.27)

pheatmap(ac.tx, 
         color = magma(100),
         breaks = breaksList,
         labels_col = paste0(colnames(ac.tx), " (", cnt.taxa$unique_gotus[match(colnames(ac.tx), cnt.taxa$phylum)], ")"),
         cutree_rows = 4,
         fontsize = 16)

dev.off()

