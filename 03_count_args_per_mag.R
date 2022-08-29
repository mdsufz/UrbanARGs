# This script plots a histogram of the number of ARGs per gOTU clustered MAGs

# packages
library(data.table)
library(ggplot2)
library(dplyr)

# load gOTU data, only want to analyze the MAGs that got clustered
gotu.final <- fread("data\\gotu_final.tsv", sep = "\t")

# only ANI 95
gotu.final <- gotu.final[which(gotu.final$ani == 95), ]

# load ARG data
deeparg <- fread("data\\deeparg.tsv", sep = "\t")

# only consider ARGs with 35% or higher percent identity
deeparg <- deeparg[deeparg$identity >= 35, ]

# count number of unique args per mag
mag.arg <- merge(gotu.final[, c("bin", "new_group")], deeparg[, c("mag", "arg")], by.x = "bin", by.y = "mag", all.x = T)
mag.arg.cnt <- as.data.frame(
  mag.arg %>% group_by(bin) %>% summarize(unique.args = n_distinct(arg, na.rm = T))
)

# create histogram with number of args per mag
pdf(file = paste0("figures\\arg_count_histogram_gOTUs.pdf"),
    height = 8.27, width = 11.69)

ggplot(mag.arg.cnt, aes(x = unique.args)) + 
  geom_histogram(binwidth = 1, color = "darkcyan", fill = "darkturquoise") + 
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    text = element_text(size = 20)
    ) +
  labs(x = "Number of ARGs per MAG", y = "Number of MAGs")

dev.off()

# average number of ARGs per gOTU MAGs
mean(mag.arg.cnt$unique.args)
sd(mag.arg.cnt$unique.args)
