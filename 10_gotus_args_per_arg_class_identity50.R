# This script plots a scatter plot the number of gOTU clusters containging an 
# ARG in each ARG class as function of the number of ARGs per ARG class
# with percent identity over 50

# packages
library(data.table)
library(ggplot2)
library(dplyr)
library(viridis)
library(forcats)

# load data
# load gOTU data, only want to analyze the MAGs that got clustered
gotu.final <- fread("data\\final_groups_qual.tsv", sep = "\t")
# only ANI 95
gotu.final <- gotu.final[which(gotu.final$ani == 95), ]

# load ARG data
deeparg <- fread("data\\deeparg.tsv", sep = "\t")
# only consider ARGs with 50% or higher percent identity
deeparg <- deeparg[deeparg$identity >= 50, ]

# count number of gOTUs per ARG class and count number of ARGs per ARG class
tmp <- merge(gotu.final[, c("bin", "new_group")], deeparg[, c("mag", "predicted.arg.class", "arg")],
             by.x = "bin", by.y = "mag", all.x = T)
arg.class.gotu.cnt <- as.data.frame(
  tmp %>% group_by(predicted.arg.class) %>% summarize(arg.cnt = n_distinct(arg), gotu.cnt = n_distinct(new_group))
)
# remove NA arg class
arg.class.gotu.cnt <- arg.class.gotu.cnt[!is.na(arg.class.gotu.cnt$predicted.arg.class),]

# divide number of unique args by number of gOTU clusters
arg.class.gotu.cnt$var <- arg.class.gotu.cnt$arg.cnt/length(unique(gotu.final$new_group))
# arg.class.gotu.cnt$var <- arg.class.gotu.cnt$gotu.cnt/arg.class.gotu.cnt$arg.cnt

# number of gotus with arg
tmp1 <- as.data.frame(tmp %>% group_by(predicted.arg.class, arg) %>% summarize(n.gotus = n_distinct(new_group)))
tmp1.red <- tmp1[tmp1$predicted.arg.class %in% arg.class.gotu.cnt$predicted.arg.class[which(arg.class.gotu.cnt$gotu.cnt > 50)],]

arg.cnt <- as.data.frame(tmp1 %>% group_by(predicted.arg.class) %>% summarise(n.args = n_distinct(arg)))
tmp1.red$arg.class.cnt <- paste0(arg.cnt$predicted.arg.class[match(tmp1.red$predicted.arg.class, arg.cnt$predicted.arg.class)], 
                                 " (", arg.cnt$n.args[match(tmp1.red$predicted.arg.class, arg.cnt$predicted.arg.class)], ")")

pdf(file = paste0("figures\\identity50_bar_gotus_per_arg.pdf"),
    width = 8.27, height = 11.69)
ggplot(tmp1.red, aes(x = fct_rev(fct_infreq(as.factor(arg.class.cnt), n.gotus)), y = n.gotus, fill = predicted.arg.class)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  labs(x = NULL, y = "Number of gOTUs with ARG") +
  theme_bw() +
  theme(text = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  coord_flip() 
dev.off()

# identify highly prevalent ARGs
tmp2 <- tmp1.red
tmp2 <- merge(tmp2, deeparg[, c("arg", "best.hit")], by = "arg", all.x = T)
tmp2 <- tmp2[order(tmp2$n.gotus, decreasing = T),]
unique(tmp2[tmp2$n.gotus >= 50, c("arg", "predicted.arg.class", "n.gotus")])
