# This script plots a bar plot of the number of gOTU species with at least one
# ARG in each ARG class

# packages
library(data.table)
library(ggplot2)
library(dplyr)
library(viridis)

# load gOTU data, only want to analyze the MAGs that got clustered
gotu.final <- fread("data\\gotu_final.tsv", sep = "\t")

# only ANI 95
gotu.final <- gotu.final[which(gotu.final$ani == 95), ]

# load ARG data
deeparg <- fread("data\\deeparg.tsv", sep = "\t")

# only consider ARGs with 35% or higher percent identity
deeparg <- deeparg[deeparg$identity >= 35, ]

# count number of unique gOTU species per ARG class
gotu.arg.class <- merge(gotu.final[, c("bin", "new_group")], deeparg[, c("mag", "arg")], by.x = "bin", by.y = "mag", all.x = T)
gotu.arg.class.cnt <- as.data.frame(
  gotu.arg.class %>% group_by(predicted.arg.class) %>% summarize(unique.gotus = n_distinct(new_group, na.rm = T))
)

# add prevalence for coloring scale
gotu.arg.class.cnt$prevalence <- gotu.arg.class.cnt$unique.gotus/length(unique(gotu.final$new_group))

# create a bar chart with number of gOTU species per ARG class
pdf(file = paste0("figures\\gOTUs_per_arg_class.pdf"),
    width = 8.27, height = 11.69)

ggplot(gotu.arg.class.cnt, aes(x = reorder(predicted.arg.class, unique.gotus), y = unique.gotus, fill = cut(prevalence, 100))) + 
  geom_col() +
  scale_fill_viridis(begin = 0, end = 1, discrete = T, option = "magma")+
  scale_y_continuous(limits = c(0, length(unique(gotu.final$new_group)))) +
  labs(x = NULL, y = "Number of gOTU clusters") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

dev.off()
