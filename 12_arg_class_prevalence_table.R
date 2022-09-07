# This script gives out a table with the ARG class prevalence per gOTU cluster

# packages
library(data.table)
library(dplyr)
library(stringr)

# load data
# load gOTU data, only want to analyze the MAGs that got clustered
gotu.final <- fread("data\\final_groups_qual.tsv", sep = "\t")
# only ANI 95
gotu.final <- gotu.final[which(gotu.final$ani == 95), ]

# load ARG data
deeparg <- fread("data\\deeparg.tsv", sep = "\t")
# only consider ARGs with 35% or higher percent identity
deeparg <- deeparg[deeparg$identity >= 35, ]

# merge data
tmp <- merge(gotu.final, deeparg, by.x = "bin", by.y = "mag", all.x = T)

# summarize data
arg.class.gotu.cnt <- 
  merge(
    as.data.frame(tmp %>% 
      group_by(new_group, predicted.arg.class) %>% 
      summarize(
        gotu.arg.class.count = n_distinct(bin, na.rm = T))
      ),
    as.data.frame(tmp %>% 
      group_by(new_group) %>% 
      summarize(
        gotu.count = n_distinct(bin, na.rm = T))
      ),
    by = "new_group"
    
  )
arg.class.gotu.cnt$prevalence <- arg.class.gotu.cnt$gotu.arg.class.count/arg.class.gotu.cnt$gotu.count
arg.class.gotu.cnt <- arg.class.gotu.cnt[!is.na(arg.class.gotu.cnt$predicted.arg.class),]

