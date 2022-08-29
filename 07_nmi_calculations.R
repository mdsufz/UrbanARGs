# This script calculates NMI and p-value of all taxa-ARG class cmobinations

# packages
library(data.table)
library(dplyr)
library(stringr)
library(aricode) # NMI

# load data
# load gOTU data, only want to analyze the MAGs that got clustered
gotu.final <- fread("data\\gotu_final.tsv", sep = "\t")
# only ANI 95
gotu.final <- gotu.final[which(gotu.final$ani == 95), ]

# load ARG data
deeparg <- fread("data\\deeparg.tsv", sep = "\t")
# only consider ARGs with 35% or higher percent identity
deeparg <- deeparg[deeparg$identity >= 35, ]

# load best bins
gotu.best <- fread("data\\gotu_best.tsv", sep = "\t")
gotu.best <- gotu.best[which(gotu.best$ani == 95), ]

# make sure taxonomy is consistent per cluster
taxonomy <- merge(gotu.final[, c("bin", "new_group")], gotu.best[, c("new_group", "taxonomy")], by = "new_group", all.x = T)
taxonomy <- as.data.frame(taxonomy)

# extract taxonomic levels
# remove "s__" etc. (e.g. "p__Firmicutes" to "Firmicutes")
taxonomy$domain <- str_extract(taxonomy$taxonomy, "(?<=d__).*?(?=;)")
taxonomy$phylum <- str_extract(taxonomy$taxonomy, "(?<=p__).*?(?=;)")
taxonomy$class <- str_extract(taxonomy$taxonomy, "(?<=c__).*?(?=;)")
taxonomy$order <- str_extract(taxonomy$taxonomy, "(?<=o__).*?(?=;)")
taxonomy$family <- str_extract(taxonomy$taxonomy, "(?<=f__).*?(?=;)")
taxonomy$genus <- str_extract(taxonomy$taxonomy, "(?<=g__).*?(?=;)")
taxonomy$species <- str_extract(taxonomy$taxonomy, "(?<=s__).*")

# remove alphabetic suffixes introduced by GTDB
# "_A" etc. (e.g. "Firmicutes_A" to "Firmicutes"), sometimes 2 letters, e.g. "_AB"
taxonomy$phylum <- gsub("_([[:alpha:]]{1,2})", "", taxonomy$phylum)
taxonomy$class <- gsub("_([[:alpha:]]{1,2})", "", taxonomy$class)
taxonomy$order <- gsub("_([[:alpha:]]{1,2})", "", taxonomy$order)
taxonomy$family <- gsub("_([[:alpha:]]{1,2})", "", taxonomy$family)
taxonomy$genus <- gsub("_([[:alpha:]]{1,2})", "", taxonomy$genus)
taxonomy$species <- gsub("_([[:alpha:]]{1,2})", "", taxonomy$species)

# replace empty with NA
taxonomy$phylum[taxonomy$phylum == ""] <- NA
taxonomy$class[taxonomy$class == ""] <- NA
taxonomy$order[taxonomy$order == ""] <- NA
taxonomy$family[taxonomy$family == ""] <- NA
taxonomy$genus[taxonomy$genus == ""] <- NA
taxonomy$species[taxonomy$species == ""] <- NA


# calculate NMI and chi-squared p-value per taxonomic category per ARG class
arg.classes <- unique(deeparg$predicted.arg.class)

# set up data frame to collect results
nmi.results <- rbind(
  setNames(do.call("rbind", replicate(length(arg.classes), as.data.frame(unique(cbind("domain", taxonomy$domain))), simplify = F)), c("tax.level", "taxon")),
  setNames(do.call("rbind", replicate(length(arg.classes), as.data.frame(unique(cbind("phylum", taxonomy$phylum))), simplify = F)), c("tax.level", "taxon")),
  setNames(do.call("rbind", replicate(length(arg.classes), as.data.frame(unique(cbind("class", taxonomy$class))), simplify = F)), c("tax.level", "taxon")),
  setNames(do.call("rbind", replicate(length(arg.classes), as.data.frame(unique(cbind("order", taxonomy$order))), simplify = F)), c("tax.level", "taxon")),
  setNames(do.call("rbind", replicate(length(arg.classes), as.data.frame(unique(cbind("family", taxonomy$family))), simplify = F)), c("tax.level", "taxon")),
  setNames(do.call("rbind", replicate(length(arg.classes), as.data.frame(unique(cbind("genus", taxonomy$genus))), simplify = F)), c("tax.level", "taxon")),
  make.row.names = F)

# make sure that there are unique combinations of taxons and arg classes
nmi.results$arg.class <- c(
  rep(arg.classes, each = length(unique(taxonomy$domain))),
  rep(arg.classes, each = length(unique(taxonomy$phylum))),
  rep(arg.classes, each = length(unique(taxonomy$class))),
  rep(arg.classes, each = length(unique(taxonomy$order))),
  rep(arg.classes, each = length(unique(taxonomy$family))),
  rep(arg.classes, each = length(unique(taxonomy$genus))))
# NA values for calculated numbers
nmi.results$nmi <- NA
nmi.results$p.value <- NA
nmi.results$n <- NA

# combine all data
myData <- merge(gotu.final, deeparg, by.x = "bin", by.y = "mag", all.x = T)
myData <- merge(myData, taxonomy, by = "bin", all.x = T)
myData <- as.data.frame(myData)

# calculate NMI and p-value
for (idx in 1:nrow(nmi.results)) {
  # extract relevant data
  extr.data <- myData[which((myData[, colnames(myData) == nmi.results$tax.level[idx]] == nmi.results$taxon[idx]) & 
                              (myData$predicted.arg.class == nmi.results$arg.class[idx])), 
                      c("arg", colnames(myData)[which(colnames(myData) == nmi.results$tax.level[idx]) + 1])]
  
  # remove NA data
  extr.data <- extr.data[complete.cases(extr.data),]
  
  # if the both variables have different level values more than two
  if (nrow(extr.data) >= 1 & length(unique(extr.data[, 1])) > 1 & length(unique(extr.data[, 2]))> 1 ) {
    # compute NMI for all ARG genes in the ARG class in all species found in the taxa
    nmi.results$nmi[idx] <- NMI(extr.data[, 1], extr.data[, 2])
    nmi.results$n[idx] <- nrow(extr.data)
    nmi.results$p.value[idx] <- chisq.test(extr.data[, 1], extr.data[, 2])$p.value
  }
}

# extract the significant results
nmi.report <- nmi.results[which((nmi.results$nmi > 0.5) & (nmi.results$p.value < 0.05)), ]

# save tables
fwrite(nmi.results, "tables\\nmi_results.csv", sep = ",", row.names = F)
fwrite(nmi.report, "tables\\nmi_significant.csv", sep = ",", row.names = F)

# plot the significant NMI ARGs
for (i in 1:nrow(nmi.report)) {
  # extract data
  extr.data <- myData[which((myData[which(colnames(myData) == nmi.report$tax.level[i])] ==  nmi.report$taxon[i]) & 
                 (myData$predicted.arg.class == nmi.report$arg.class[i])), ]
  
  # dplyr needs a specific column name, can't use a variable
  extr.data$taxons <- extr.data[, which(colnames(extr.data) == nmi.report$tax.level[i]) + 1]
  
  # data to plot
  data.plot <- as.data.frame(
    extr.data %>% group_by(taxons, arg) %>% summarize(n = n_distinct(bin))
  )
  
  # plot
  pdf(file = paste0("figures\\nmi_sign_", nmi.report$taxon[i], "_", nmi.report$arg.class[i], ".pdf"),
      height = 8.27, width = 11.69)
  print(ggplot(data.plot, aes(x = taxons, y = n, fill = arg)) +
    geom_col(position="stack") +
    theme_bw() + 
    theme(text = element_text(size = 18),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    labs(x = colnames(extr.data[which(colnames(extr.data) == nmi.report$tax.level[i]) + 1]), y = "Number of MAGs") +
    coord_flip() +
    guides(fill = guide_legend(title = "ARG")))
  dev.off()
}



