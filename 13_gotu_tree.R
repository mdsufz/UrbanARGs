library(data.table)
library(stringi)
library(ggplot2)
library(ggtree)
library(dplyr)
library(viridis)

## create tree, 
# load best bins to assign taxonomy to clusters
gotu.best <- fread("data\\bestbins_metadata.tsv", sep = "\t")
gotu.best <- gotu.best[gotu.best$ani == 95, ]

# get representatives (last non-emtpy identifier)
gotu.best.rep <- stri_extract_last(gotu.best$taxonomy, regex="(?<=(^|;))([dpcofgs]{1}_{2}[a-zA-Z0-9_\\- ]+)(?=(;*))")
# fwrite(as.list(gotu.best.rep), "trees\\gotu.representatives.txt", sep = ",")

# ids that are not found in phyloT and the replacements
# g__GCA-016699365
gotu.best.rep[gotu.best.rep == "g__GCA-016699365"] <- "f__UBA9973"
# g__JAGOQQ01
gotu.best.rep[gotu.best.rep == "g__JAGOQQ01"] <- "f__UBA9973"
# g__JACRIG01
gotu.best.rep[gotu.best.rep == "g__JACRIG01"] <- "f__Peribacteraceae"
# g__UBA11153
gotu.best.rep[gotu.best.rep == "g__UBA11153"] <- "f__Coleofasciculaceae"
# f__DSSB01
gotu.best.rep[gotu.best.rep == "f__DSSB01"] <- "o__DSSB01"
# g__JAGORN01
gotu.best.rep[gotu.best.rep == "g__JAGORN01"] <- "f__CAIOMD01"
# g__UBA4116
gotu.best.rep[gotu.best.rep == "g__UBA4116"] <- "f__Tannerellaceae"

# fwrite(as.list(gotu.best.rep), "trees\\gotu.representatives.txt", sep = ",")
# create tree with phyloT
# load tree
gotu.tree <- treeio::read.newick("trees\\tree_gotu_newick.txt")

# basic tree
# ggtree(gotu.tree, layout="circular", branch.length='none')

# plot prevalences
# load gOTU data, only want to analyze the MAGs that got clustered
gotu.final <- fread("data\\final_groups_qual.tsv", sep = "\t")
# only ANI 95
gotu.final <- gotu.final[which(gotu.final$ani == 95), ]

# load ARG data
deeparg <- fread("data\\deeparg.tsv", sep = "\t")
# only consider ARGs with 35% or higher percent identity
deeparg <- deeparg[deeparg$identity >= 35, ]

# load best bins to assign taxonomy to clusters
gotu.best <- fread("data\\bestbins_metadata.tsv", sep = "\t")
gotu.best <- gotu.best[gotu.best$ani == 95, ]
# put all data together in a table
arg.class.tax <- merge(gotu.final[, c("bin", "new_group")], deeparg[, c("mag", "arg", "predicted.arg.class")], 
                       by.x = "bin", by.y = "mag", all.x = T)
# add tree node label
gotu.best$tree.id <- gotu.best.rep
arg.class.tax$leaf <- gsub(" ", "_", gotu.best$tree.id[match(arg.class.tax$new_group, best.bins$new_group)])


# summarize data
leaf.prevalence <- 
  merge(
    as.data.frame(arg.class.tax %>% 
                    group_by(leaf, predicted.arg.class) %>% 
                    summarize(
                      gotu.arg.class.count = n_distinct(bin))
    ),
    as.data.frame(arg.class.tax %>% 
                    group_by(leaf) %>% 
                    summarize(
                      mag.count = n_distinct(bin))
    ),
    by = "leaf"
    
  )
# remove arg class NA
leaf.prevalence <- leaf.prevalence[!is.na(leaf.prevalence$predicted.arg.class),]
# calculate prevalence
leaf.prevalence$prevalence <- leaf.prevalence$gotu.arg.class.count/leaf.prevalence$mag.count

# plot prevalence heatmap onto tree
leaf.matrix <- reshape2::dcast(leaf.prevalence, leaf~predicted.arg.class, value.var = "prevalence")
leaf.matrix <- leaf.matrix[!is.na(leaf.matrix$leaf),]
rownames(leaf.matrix) <- leaf.matrix$leaf
leaf.matrix$leaf <- NULL
leaf.matrix[is.na(leaf.matrix)] <- 0
# add leaves that have 0 prevalence
gotu.best$leaf <- gsub(" ", "_", gotu.best$tree.id)
missing <- gotu.best$leaf[!(gotu.best$leaf %in% rownames(leaf.matrix))]
missing.mat <- matrix(0, nrow = length(missing), ncol = ncol(leaf.matrix))
rownames(missing.mat) <- missing
colnames(missing.mat) <- colnames(leaf.matrix)
leaf.matrix <- rbind(leaf.matrix, missing.mat)

# only include arg classes, skip absent class
# use WAP, calculated in script 05
leaf.matrix <- leaf.matrix[, rev(names(wap)[wap > 0.25])]
# use NMI results, calculated in script 07

paste0(substr(nmi.report$tax.level, 1, 1), "__", nmi.report$taxon)

d <- as_tibble(gotu.tree)
tr3 <- tidytree::as.treedata(d)

dev.off()
# plot the tree with heatmap
pdf(file = paste0("figures\\gotu_tree.pdf"),
    height = 7.87, width = 7.87)
  p.tree <- ggtree(gotu.tree, layout="fan", open.angle = 20, branch.length='none')
  
  p.tree <- p.tree +
    geom_hilight(node=which(p.tree$data$label == "p__Proteobacteria"), fill="deeppink", alpha=.25) +
    geom_hilight(node=which(p.tree$data$label == "p__Bacteroidota"), fill="darkorchid1", alpha=.25) +
    geom_hilight(node=which(p.tree$data$label == "p__Firmicutes"), fill="goldenrod1", alpha=.25) +
    geom_hilight(node=which(p.tree$data$label == "p__Firmicutes_A"), fill="goldenrod1", alpha=.25) +
    geom_hilight(node=which(p.tree$data$label == "p__Firmicutes_C"), fill="goldenrod1", alpha=.25) +
    geom_hilight(node=which(p.tree$data$label == "p__Cyanobacteria"), fill="darkgreen", alpha=.25) +
    geom_hilight(node=which(p.tree$data$label == "p__Actinobacteriota"), fill="steelblue1", alpha=.25) +
    geom_hilight(node=which(p.tree$data$label == "p__Patescibacteria"), fill="salmon4", alpha=.25)
  
  p.tree <- p.tree + 
    geom_point2(aes(subset = (node %in% which(p.tree$data$label %in% paste0(substr(nmi.report$tax.level, 1, 1), "__", nmi.report$taxon)))),
      size = 3, colour = 'black', fill = "yellow", shape = 21)
  # p.tree <- p.tree + geom_label2(aes(subset = (node %in% which(p.tree$data$label %in% paste0(substr(nmi.report$tax.level, 1, 1), "__", nmi.report$taxon))), label = label))
  
  # add heatmap
  gheatmap(p.tree, leaf.matrix, color = NULL, colnames = T, width = 0.67, legend_title = "Prevalence") +
    scale_fill_gradientn(limits = c(0, 1), colours=viridis(10,o="A")) +
    theme(legend.position = "none")
dev.off()
