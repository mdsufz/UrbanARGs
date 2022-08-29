# This script plots a scatter plot the number of gOTU clusters containging an 
# ARG in each ARG class as function of the number of ARGs per ARG class

# packages
library(data.table)
library(ggplot2)
library(dplyr)
library(viridis)
library(forcats)

# load data
# load gOTU data, only want to analyze the MAGs that got clustered
gotu.final <- fread("data\\gotu_final.tsv", sep = "\t")
# only ANI 95
gotu.final <- gotu.final[which(gotu.final$ani == 95), ]

# load ARG data
deeparg <- fread("data\\deeparg.tsv", sep = "\t")
# only consider ARGs with 35% or higher percent identity
deeparg <- deeparg[deeparg$identity >= 35, ]

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

# scatter plot
pdf(file = paste0("figures\\args_gotus_per_arg_class.pdf"),
    width = 8.27, height = 11.69)

ggplot(arg.class.gotu.cnt, aes(x = arg.cnt, y = gotu.cnt, label = predicted.arg.class)) +
  geom_point(size = 4) +
  theme_bw() +
  theme(text = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  labs(x = "Number of ARGs", y = "Number of gOTU clusters") +
  geom_text(aes(label = ifelse(gotu.cnt > 75, as.character(predicted.arg.class), '')), hjust=-0.15, vjust=0)

dev.off()

# calculate Kendall's tau and chi-squared p-value for correlation
# between all arg classes
ktau.all <- cor(arg.class.gotu.cnt$arg.cnt, arg.class.gotu.cnt$gotu.cnt, method = "kendall")
# between highly prevalent arg classes
ktau.prev <- cor(arg.class.gotu.cnt$arg.cnt[arg.class.gotu.cnt$gotu.cnt > 100], 
            arg.class.gotu.cnt$gotu.cnt[arg.class.gotu.cnt$gotu.cnt > 100], method = "kendall")
# between very highly prevalent arg classes
ktau.high <- cor(arg.class.gotu.cnt$arg.cnt[arg.class.gotu.cnt$gotu.cnt > 200], 
                 arg.class.gotu.cnt$gotu.cnt[arg.class.gotu.cnt$gotu.cnt > 200], method = "kendall")

# bar plot of the number of args dividied by gotus
# only for arg classes that appear in more than 100 gotus
# the others may also be variable, but there are too few data points to tell
pdf(file = paste0("figures\\bar_args_gotus_per_arg_class.pdf"),
    width = 8.27, height = 11.69)

# only plot higher than 100 gotus
df.red <- arg.class.gotu.cnt[which(arg.class.gotu.cnt$gotu.cnt > 100),]
ggplot(df.red, aes(x = reorder(predicted.arg.class, var), y = var)) +
  geom_col() +
  theme_bw() +
  theme(text = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  labs(x = NULL, y = "Ratio of ARGs per gOTUs") +
  coord_flip()

dev.off()

# number of gotus with arg
tmp1 <- as.data.frame(tmp %>% group_by(predicted.arg.class, arg) %>% summarize(n.gotus = n_distinct(new_group)))
tmp1.red <- tmp1[tmp1$predicted.arg.class %in% arg.class.gotu.cnt$predicted.arg.class[which(arg.class.gotu.cnt$gotu.cnt > 100)],]

arg.cnt <- as.data.frame(tmp1 %>% group_by(predicted.arg.class) %>% summarise(n.args = n_distinct(arg)))
tmp1.red$arg.class.cnt <- paste0(arg.cnt$predicted.arg.class[match(tmp1.red$predicted.arg.class, arg.cnt$predicted.arg.class)], 
                                 " (", arg.cnt$n.args[match(tmp1.red$predicted.arg.class, arg.cnt$predicted.arg.class)], ")")

pdf(file = paste0("figures\\bar_gotus_per_arg.pdf"),
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
unique(tmp2[tmp2$n.gotus >= 100, c("arg", "predicted.arg.class", "n.gotus")])
