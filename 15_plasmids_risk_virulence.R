setwd("Y:\\Home\\magnusdo\\projects\\urban_mags_args")

# packages
library(data.table)
library(ggplot2)
library(dplyr)
library(viridis)
library(pheatmap)
library(stringr)

# load data
# load gOTU data, only want to analyze the MAGs that got clustered
gotu.final <- fread("data\\final_groups_qual.tsv", sep = "\t")
# only ANI 95
gotu.final <- gotu.final[which(gotu.final$ani == 95), ]

# load best bins
gotu.best <- fread("data\\bestbins_metadata.tsv", sep = "\t")
gotu.best <- gotu.best[gotu.best$ani == 95, ]

# load ARG data
deeparg <- fread("data\\deeparg.tsv", sep = "\t")
# only consider ARGs with 35% or higher percent identity
deeparg <- deeparg[deeparg$identity >= 35, ]
deeparg <- deeparg[deeparg$mag %in% gotu.final$bin,]

# merge data
all.data <- merge(gotu.final[, c("bin", "new_group")], 
                  deeparg[, c("mag", "arg", "predicted.arg.class", "read.id")], 
                  by.x = "bin", by.y = "mag", all = T)

# load plasflow data
plasflow <- fread("data\\plasflow.tsv", sep = "\t", header = T)
otu.map <- fread("Y:\\Home\\magnusdo\\projects\\urban\\CTbins-tracking-table.tsv", sep = "\t", header = F)
plasflow$mag <- otu.map$V2[match(gsub(".tsv", ".fa", plasflow$EASRX2001224_bin.10.tsv), otu.map$V1)]
plasflow$mag <- gsub(".fa", "", plasflow$mag)
plasflow$read.type <- gsub("\\.[[:alnum:]]+", "", plasflow$label)

# add read type
all.data$read.type <- plasflow$read.type[match(all.data$read.id, plasflow$contig_name)]

# load virulence data
virulence <- fread("data\\blastx.out", header = F, sep = "\t")
colnames(virulence) <- c("qseqid", "sseqid", "pident", "length", "mismatch", 
                         "gapopen", "qstart", "qend", "sstart", "send", "evalue", 
                         "bitscore", "qcovs", "qcovhsp", "mag")
virulence$mag <- gsub(".out", ".fa", virulence$mag)
virulence$mag <- otu.map$V2[match(virulence$mag, otu.map$V1)]
virulence$mag <- gsub(".fa", "", virulence$mag)
virulence <- virulence[which(virulence$pident >= 80),]

all.data <- merge(all.data, virulence[, c("mag", "qseqid", "sseqid")], 
                  by.x = c("bin", "read.id"), 
                  by.y = c("mag", "qseqid"), all.x = T)

#
# load risk assessment from Zhang et al 2021 ()
#risk.assessment <- fread("data\\Zhang2021_SupplementaryData4.csv", sep = "\t", header = T)

#deeparg$gene.id <- str_extract(deeparg$best.hit, "[^|]+")
#deeparg$gene.id <- gsub("gi:[0-9]+:[[:alnum:]]+:", "", deeparg$gene.id)
#deeparg$gene.id <- gsub(":[[:alnum:]]*", "", deeparg$gene.id)
#gene.id <- unique(deeparg$gene.id)

#deeparg$risk.rank <- risk.assessment$Rank[match(deeparg$gene.id, risk.assessment$ARG)]

# library(treemapify)
# library(formattable)
# deeparg[(deeparg$predicted.arg.class %in% 
#           c("multidrug", "glycopeptide", "beta-lactam", "MLS", "diaminopyrimidine",
#             "bacitracin", "aminoglycoside", "peptide", "phenicol", "fluoroquinolone",
#             "sulfonamide", "tetracycline", "unclassified")) & 
#           (!is.na(deeparg$arg)) &
#           (!is.na(deeparg$read.type)),] %>% 
#   group_by(predicted.arg.class, read.type) %>%
#   summarize(n = n_distinct(arg)) %>%
#   mutate(freq =formattable::percent(n / sum(n), digits = 0)) %>%
#   ggplot(., aes(area = n, fill = read.type, subgroup = predicted.arg.class, label = freq)) +
#   scale_fill_manual(values=c("#7e4e90ff", "#de7065ff", "#f9b641ff")) +
#   geom_treemap() +
#   geom_treemap_subgroup_border(colour = "white", size = 5) +
#   geom_treemap_subgroup_text(place = "bottom", grow = TRUE,
#                              alpha = 1, colour = "black",
#                              fontface = "italic") +
#   geom_treemap_text(colour = "white",
#                     place = "topright",
#                     size = 16)
# 
# deeparg$risk.rank[is.na(deeparg$risk.rank)] <- "Unassessed"
# deeparg$risk.rank <- paste0("Risk rank: ", deeparg$risk.rank)
# deeparg$num.virulence[is.na(deeparg$num.virulence)] <- 0
# #deeparg$num.virulence[deeparg$num.virulence == 0] <- NA
# #deeparg <- mutate(deeparg, log.num.vir = num.virulence + 0.01)



# How many ARGs are on chromosomes and plasmids?
num.arg.read.type <- unique(all.data[which(!is.na(all.data$arg)), c("bin", "read.id", "arg", "read.type")]) %>% 
  group_by(read.type) %>% 
  count(read.type)
num.arg.read.type$perc <- round(100*(num.arg.read.type$n/sum(num.arg.read.type$n)),0)

# pdf(file = "figures\\num.arg.read.type.pdf", width = 11.7, height = 8.3)
# ggplot(subset(num.arg.read.type, !is.na(read.type)), aes(x = "", y = n, fill = read.type)) + 
#   geom_bar(width = 1, stat = "identity", color = "black") + 
#   coord_polar("y", start = 0) + 
#   scale_fill_manual(values = c("#0073C2FF", "#CD534CFF", "#EFC000FF")) +
#   theme(axis.text.x = element_blank(), 
#         panel.background = element_rect(fill = 'white')
#         ) +
#   labs(x = NULL, y = NULL)
# dev.off()
# 
# # How many virulent factors align to each gOTU species?
# 
# # at all
# pdf(file = "figures\\virulence_gotus.pdf", width = 6, height = 4.5)
# all.data[, c("new_group", "sseqid")] %>% group_by(new_group) %>% summarise(n = n_distinct(sseqid, na.rm = T)) %>%
#   ggplot(., aes(x = n)) + 
#   geom_histogram(bins = 30, color = "black", fill = "#FF6666") +
#   # geom_vline(xintercept = 100, linetype="longdash", color = "darkgrey") +
#   theme_bw() +
#   theme(text = element_text(size = 16)) +
#   labs(x = "Number of aligned virulence factors", y = "Number of gOTUs")
# dev.off()
# 
# # per read type
# pdf(file = "figures\\virulence_gotus_read_type.pdf", width = 11.7, height = 4.5)
# all.data[!is.na(all.data$read.type), c("new_group", "sseqid", "read.type")] %>% group_by(new_group, read.type) %>% summarise(n = n_distinct(sseqid, na.rm = T)) %>%
#   ggplot(., aes(x = n)) + 
#   geom_histogram(bins = 30, color = "black", aes(fill = read.type)) +
#   theme_bw() +
#   theme(text = element_text(size = 16), legend.position = "none") +
#   labs(x = "Number of aligned virulence factors", y = "Number of gOTUs") +
#   facet_grid(. ~ read.type) +
#   scale_fill_manual(values = c("#0073C2FF", "#CD534CFF", "#EFC000FF"))
# dev.off()

# how many reads align at all, do mean and sd
num.vir.gotu <- all.data[, c("new_group", "sseqid")] %>% 
  group_by(new_group) %>% 
  summarise(n = n_distinct(sseqid, na.rm = T))
# mean(num.vir.gotu$n)
# sd(num.vir.gotu$n)
# median(num.vir.gotu$n)
# min(num.vir.gotu$n)
# max(num.vir.gotu$n)
# 
# # number of virulence factors per gOTU with args in each arg class
# # widespread and ubiquitous
# pdf(file = "figures\\virulence_arg_class_wu.pdf", width = 6, height = 6)
# all.data[!is.na(all.data$arg) & (all.data$predicted.arg.class %in% c("multidrug", "glycopeptide", "MLS", "bacitracin", "beta-lactam", "unclassified")), 
#          c("new_group", "sseqid", "arg", "predicted.arg.class")] %>%
#   group_by(predicted.arg.class, new_group) %>% 
#   summarise(n = n_distinct(sseqid, na.rm = T)) %>%
#   #filter(any(n > 100)) %>%
#   ggplot(., aes(x = n, y = reorder(predicted.arg.class, n, FUN = max))) +
#     geom_jitter(aes(color = predicted.arg.class), alpha = 0.5, size = 2) +
#     geom_boxplot(outlier.shape = NA, alpha = 0) +
#     scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 10,  100, 200, 400, 600)) +
#     theme_bw() +
#     theme(legend.position = "none", text = element_text(size = 16)) +
#     labs(x = "Number of virulence factors per gOTU", y = NULL)
# dev.off()
# # all
# pdf(file = "figures\\virulence_arg_class_all.pdf", width = 11.7, height = 8.3)
# all.data[!is.na(all.data$arg), c("new_group", "sseqid", "arg", "predicted.arg.class")] %>%
#   group_by(predicted.arg.class, new_group) %>% 
#   summarise(n = n_distinct(sseqid, na.rm = T)) %>%
#   #filter(any(n > 100)) %>%
#   ggplot(., aes(x = n, y = reorder(predicted.arg.class, n, FUN = max))) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(aes(color = predicted.arg.class), alpha = 0.75, size = 2) +
#   theme(legend.position = "none", text = element_text(size = 16)) +
#   theme_bw() +
#   labs(x = "Number of virulence factors per gOTU", y = NULL)
# dev.off()

# pdf(file = "figures\\risk_rank_read_type_virulence.pdf", width = 11.7, height = 8.3)
# ggplot(deeparg[(deeparg$predicted.arg.class %in% 
#                   c("multidrug", "glycopeptide", "beta-lactam", "MLS", "bacitracin", "unclassified")) & 
#                  (!is.na(deeparg$arg)) & 
#                  (!is.na(deeparg$read.type)),], 
#        aes(x = predicted.arg.class, y = num.virulence, fill = predicted.arg.class)) + 
#   #geom_violin(alpha=0.5) + 
#   geom_jitter(shape = 21, size=3, alpha=1, width = .2) +
#   #scale_y_continuous(trans="log10") +
#   coord_flip() +
#   facet_grid(read.type ~ risk.rank) +
#   labs(y = "Number of virulence factors found on same scaffold as ARG", x = NULL) +
#   theme(legend.position = "none", text = element_text(size = 16))
# dev.off()
# 
# 
# mag.virulence <- as.data.frame(
#   virulence %>% group_by(mag) %>% summarize(num.virulence = n_distinct(sseqid))
# )
# 
# # histogram of number of virulence factors aligned to MAG
# ggplot(mag.virulence, aes(x = num.virulence)) + 
#   geom_histogram(binwidth = 20) + 
#   geom_vline(xintercept = 150) +
#   geom_vline(xintercept = 550)
# 
# # Classify MAGs with more than 150 aligned virulence factors as "potentially virulent"
# deeparg$mag.virulence <- "Low"
# deeparg$mag.virulence[which(deeparg$mag %in% mag.virulence$mag[mag.virulence$num.virulence > 150])] <- "Medium"
# deeparg$mag.virulence[which(deeparg$mag %in% mag.virulence$mag[mag.virulence$num.virulence > 550])] <- "High"
# 
# deeparg[!is.na(deeparg$read.type),] %>% group_by(mag.virulence, read.type) %>% summarise(n = n_distinct(arg)) %>%
#   mutate(freq =formattable::percent(n / sum(n), digits = 0)) %>%
#   ggplot(., aes(x = mag.virulence, y = n, fill = read.type)) +
#   geom_bar(position="fill", stat="identity")
# 
# 
# ggplot(deeparg[(deeparg$predicted.arg.class %in% 
#                   c("multidrug", "glycopeptide", "beta-lactam", "MLS", "bacitracin", "unclassified")) & 
#                  (!is.na(deeparg$arg)) & 
#                  (!is.na(deeparg$read.type)),], 
#        aes(x = predicted.arg.class, y = num.virulence, fill = predicted.arg.class)) + 
#   #geom_violin(alpha=0.5) + 
#   geom_jitter(shape = 21, size=3, alpha=1, width = .2) +
#   #scale_y_continuous(trans="log10") +
#   coord_flip() +
#   facet_grid(read.type ~ mag.virulence) +
#   labs(y = "Number of virulence factors found on same scaffold as ARG", x = NULL) +
#   theme(legend.position = "none", text = element_text(size = 16))
# 
# pdf(file = "figures\\virulence_multidrug.pdf", width = 8.3, height = 11.7)
# all.data[!is.na(all.data$arg) & (all.data$predicted.arg.class == "multidrug"), 
#          c("new_group", "sseqid", "arg", "predicted.arg.class")] %>%
#   group_by(arg, new_group) %>% 
#   summarise(n = n_distinct(sseqid, na.rm = T)) %>%
#   #filter(any(n > 100)) %>%
#   ggplot(., aes(x = n, y = reorder(arg, n, FUN = max))) +
#   geom_jitter(aes(color = arg), alpha = 0.5, size = 2) +
#   geom_boxplot(outlier.shape = NA, alpha = 0) +
#   scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 10,  100, 200, 400, 600)) +
#   theme_bw() +
#   theme(legend.position = "none", text = element_text(size = 16)) +
#   labs(x = "Number of virulence factors per gOTU", y = NULL)
# dev.off()
# 
# pdf(file = "figures\\virulence_multidrug_100.pdf", width = 8.3, height = 11.7)
# tmp <- all.data[!is.na(all.data$arg) & (all.data$predicted.arg.class == "multidrug"), 
#          c("new_group", "sseqid", "arg", "predicted.arg.class")] %>%
#   group_by(arg, new_group) %>% 
#   summarise(n = n_distinct(sseqid, na.rm = T)) 
# 
#   ggplot(tmp[tmp$arg %in% tmp$arg[tmp$n >= 100],], aes(x = n, y = reorder(arg, n, FUN = max))) +
#   geom_jitter(aes(color = arg), alpha = 0.5, size = 2) +
#   geom_boxplot(outlier.shape = NA, alpha = 0) +
#   scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 10,  100, 200, 400, 600)) +
#   theme_bw() +
#   theme(legend.position = "none", text = element_text(size = 16)) +
#   labs(x = "Number of virulence factors per gOTU", y = NULL)
# dev.off()

tmp.all <- all.data[!is.na(all.data$arg), c("new_group", "sseqid", "arg", "predicted.arg.class")] %>%
  group_by(predicted.arg.class, new_group) %>% 
  summarise(n = n_distinct(sseqid, na.rm = T))

pvals <- matrix(NA, nrow = length(unique(tmp.all$predicted.arg.class)), ncol = length(unique(tmp.all$predicted.arg.class)))
for (i in 1:length(unique(tmp.all$predicted.arg.class))) {
  for (j in 1:length(unique(tmp.all$predicted.arg.class))) {
    x <- tmp.all$n[which(tmp.all$predicted.arg.class == unique(tmp.all$predicted.arg.class)[i])]
    y <- tmp.all$n[which(tmp.all$predicted.arg.class == unique(tmp.all$predicted.arg.class)[j])]
    if ((length(x) >= 3) & (length(y) >= 3)) {
      pvals[i, j] <- wilcox.test(x, y)$p.value
      pvals[j, i] <- pvals[i, j]
    }
  }
}
rownames(pvals) <- unique(tmp.all$predicted.arg.class)
colnames(pvals) <- unique(tmp.all$predicted.arg.class)

pvals[rownames(pvals) %in% c("multidrug", "glycopeptide", "MLS", "bacitracin", "beta-lactam", "unclassified"), 
      colnames(pvals) %in% c("multidrug", "glycopeptide", "MLS", "bacitracin", "beta-lactam", "unclassified")]


f3a <- ggplot(subset(num.arg.read.type, !is.na(read.type)), aes(x = "", y = n, fill = read.type)) + 
  geom_bar(width = 1, stat = "identity", color = "black") + 
  coord_polar("y", start = 0) + 
  scale_fill_manual(values = c("#0073C2FF", "#CD534CFF", "#EFC000FF")) +
  theme(axis.text.x = element_blank(), 
        panel.background = element_rect(fill = 'white',)) +
  labs(x = NULL, y = NULL) +
  guides(fill="none")

f3b <- all.data[!is.na(all.data$read.type), c("new_group", "sseqid", "read.type")] %>% 
  group_by(new_group, read.type) %>% 
  summarise(n = n_distinct(sseqid, na.rm = T)) %>%
  ggplot(., aes(x = n)) + 
  geom_histogram(bins = 30, color = "black", aes(fill = read.type)) +
  theme_bw() +
  theme(text = element_text(size = 12), legend.position = "none") +
  labs(x = "Number of aligned virulence factors", y = "Number of gOTUs") +
  facet_grid(. ~ read.type) +
  scale_fill_manual(values = c("#0073C2FF", "#CD534CFF", "#EFC000FF"))


all.data[!is.na(all.data$read.type), c("new_group", "sseqid", "read.type")] %>% 
  group_by(new_group, read.type) %>% 
  summarise(n = n_distinct(sseqid, na.rm = T)) %>%
  group_by(read.type) %>%
  summarise(mean = mean(n, na.rm = T), 
                  sd = sd(n, na.rm = T), 
                  median = median(n, na.rm = T), 
                  min = min(n, na.rm = T), 
                  max = max(n, na.rm = T))

gotu.virulence <- as.data.frame(all.data[!is.na(all.data$arg), 
                                         c("new_group", "sseqid")] %>%
                                  group_by(new_group) %>% 
                                  summarise(n.virulence = n_distinct(sseqid, na.rm = T)))

gotu.virulence.argclass <- as.data.frame(all.data[!is.na(all.data$arg), 
                                         c("new_group", "predicted.arg.class", "sseqid")] %>%
                                  group_by(new_group, predicted.arg.class) %>% 
                                  summarise(n.virulence.argclass = n_distinct(sseqid, na.rm = T)))

# summarize data
gotu.prevalence <- 
  merge(
    as.data.frame(all.data %>% 
                    group_by(new_group, predicted.arg.class) %>% 
                    summarize(
                      gotu.arg.class.count = n_distinct(bin, na.rm = T))
    ),
    as.data.frame(all.data %>% 
                    group_by(new_group) %>% 
                    summarize(
                      gotu.count = n_distinct(bin, na.rm = T))
    ),
    by = "new_group"
    
  )
gotu.prevalence$prevalence <- gotu.prevalence$gotu.arg.class.count/gotu.prevalence$gotu.count
gotu.prevalence <- gotu.prevalence[!is.na(gotu.prevalence$predicted.arg.class),]

gotu.both <- merge(gotu.virulence.argclass, gotu.prevalence, by = c("new_group", "predicted.arg.class"))

gotu.both$taxonomy <- gotu.best$taxonomy[match(gotu.both$new_group, gotu.best$new_group)]
gotu.both$family <- str_extract(gotu.both$taxonomy, "(?<=f__).*?(?=;)")

# gotu.both$ARG.class <- "other"
# gotu.both$ARG.class[(gotu.both$n.virulence.argclass >= 100) & (gotu.both$prevalence >= 0.75) & (gotu.both$predicted.arg.class == "multidrug")] <- "multidrug"
# gotu.both$ARG.class[(gotu.both$n.virulence.argclass >= 100) & (gotu.both$prevalence >= 0.75) & (gotu.both$predicted.arg.class == "beta-lactam")] <- "beta-lactam"
# gotu.both$ARG.class <- as.factor(gotu.both$ARG.class)
gotu.both$ARG.class <- "other"
gotu.both$ARG.class[(gotu.both$predicted.arg.class == "multidrug")] <- "multidrug"
gotu.both$ARG.class[(gotu.both$predicted.arg.class == "beta-lactam")] <- "beta-lactam"
gotu.both$ARG.class[(gotu.both$predicted.arg.class == "glycopeptide")] <- "glycopeptide"
gotu.both$ARG.class[(gotu.both$predicted.arg.class == "MLS")] <- "MLS"
gotu.both$ARG.class[(gotu.both$predicted.arg.class == "unclassified")] <- "unclassified"
gotu.both$ARG.class[(gotu.both$predicted.arg.class == "aminoglycoside")] <- "aminoglycoside"
gotu.both$ARG.class[(gotu.both$predicted.arg.class == "tetracycline")] <- "tetracycline"
gotu.both$ARG.class <- as.factor(gotu.both$ARG.class)

# gotu.both$family.sign <- "other"
# gotu.both$family.sign[(gotu.both$n.virulence.argclass >= 100) & (gotu.both$prevalence >= 0.75) & (gotu.both$family == "Enterobacteriaceae")] <- "Enterobacteriaceae"
# gotu.both$family.sign[(gotu.both$n.virulence.argclass >= 100) & (gotu.both$prevalence >= 0.75) & (gotu.both$family == "Rhizobiaceae_A")] <- "Rhizobiaceae_A"
# gotu.both$family.sign <- as.factor(gotu.both$family.sign)
gotu.both$family.sign <- "other"
gotu.both$family.sign[(gotu.both$family == "Enterobacteriaceae")] <- "Enterobacteriaceae"
gotu.both$family.sign[(gotu.both$family == "Rhizobiaceae_A")] <- "Rhizobiaceae_A"
gotu.both$family.sign <- as.factor(gotu.both$family.sign)

f3c <- 
  # gotu.both %>%
  # group_by(predicted.arg.class, new_group) %>% 
  # summarise(n = n_distinct(sseqid, na.rm = T)) %>%
  #filter(any(n > 100)) %>%
  ggplot(gotu.both[gotu.both$predicted.arg.class %in% gotu.both$predicted.arg.class[gotu.both$n.virulence.argclass >= 100],], 
  #ggplot(gotu.both, 
         aes(x = n.virulence.argclass, 
             y = reorder(predicted.arg.class, n.virulence.argclass, FUN = max))) +
  geom_jitter(aes(fill = predicted.arg.class), 
              alpha = 0.7, 
              size = 2.5,
              shape = 21) +
  #geom_boxplot(outlier.shape = NA, 
  #             alpha = 0) +
 # scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 10,  100, 400)) +
  scale_fill_manual(values = c("#f8766dff", "#be9c00ff", "#00be70ff", "#00c1abff", "#00bbdaff", "#f962ddff", "#ff65acff")) +
  theme_bw() +
  theme(legend.position = "none", 
        text = element_text(size = 12)) +
  labs(x = "Number of unique virulence factors \n co-occurring on ARG contigs per gOTU", 
       y = NULL)


fs7 <- 
  # gotu.both %>%
  # group_by(predicted.arg.class, new_group) %>% 
  # summarise(n = n_distinct(sseqid, na.rm = T)) %>%
  #filter(any(n > 100)) %>%
  ggplot(gotu.both, 
         #ggplot(gotu.both, 
         aes(x = n.virulence.argclass, 
             y = reorder(predicted.arg.class, n.virulence.argclass, FUN = max))) +
  geom_jitter(aes(fill = predicted.arg.class), 
              alpha = 0.7, 
              size = 2.5,
              shape = 21) +
  #geom_boxplot(outlier.shape = NA, 
  #             alpha = 0) +
  # scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 10,  100, 400)) +
  #scale_fill_manual(values = c("#f8766dff", "#be9c00ff", "#00be70ff", "#00c1abff", "#00bbdaff", "#f962ddff", "#ff65acff")) +
  theme_bw() +
  theme(legend.position = "none", 
        text = element_text(size = 12)) +
  labs(x = "Number of unique virulence factors \n co-occurring on ARG contigs per gOTU", 
       y = NULL)
# tmp <- all.data[!is.na(all.data$arg) & (all.data$predicted.arg.class == "multidrug"), 
#                 c("new_group", "sseqid", "arg", "predicted.arg.class")] %>%
#   group_by(arg, new_group) %>% 
#   summarise(n = n_distinct(sseqid, na.rm = T)) 
# 
# addline_format <- function(x,...){
#   gsub('\\s','\n',x)
# }

 
  
tmp <- all.data[all.data$predicted.arg.class %in% c("multidrug", "beta-lactam", "glycopeptide", "unclassified", "tetracycline", "MLS", "aminoglycoside")] %>% 
  group_by(new_group, predicted.arg.class, arg) %>%
  summarize(n.virulence.arg = n_distinct(sseqid)) 

args <- unique(tmp$arg[tmp$arg %in% tmp$arg[tmp$n.virulence.arg >= 100]])
pvals <- matrix(NA, nrow = length(args), ncol = length(args))
for (i in 1:length(args)) {
  for (j in 1:length(args)) {
    x <- tmp$n.virulence.arg[which(tmp$arg == args[i])]
    y <- tmp$n.virulence.arg[which(tmp$arg == args[j])]
    if ((length(x) >= 3) & (length(y) >= 3)) {
      pvals[i, j] <- wilcox.test(x, y)$p.value
      pvals[j, i] <- pvals[i, j]
    }
  }
}
rownames(pvals) <- unique(args)
colnames(pvals) <- unique(args)


  # filter(any(n.virulence.arg > 100)) %>%
f3d <-  
  ggplot(tmp[tmp$arg %in% tmp$arg[tmp$n.virulence.arg >= 100],], 
         aes(x = n.virulence.arg, 
             y = reorder(arg, n.virulence.arg, FUN = max))) +
  geom_jitter(aes(fill = predicted.arg.class), 
              alpha = 0.7, 
              size = 2,
              shape = 21,
              width = 0.2) +
 # geom_boxplot(outlier.shape = NA, 
 #              alpha = 0) +
  #scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 10, 100, 200)) +
  scale_y_discrete(labels= c("smeF", "mdtG","sdiA", "multidrug \n ABC \n transporter", "vanD", "mdtM", "vanRI", "macA", "oqxB", "mdtK", "mdtC", "penA", "ompF", "ompR", "vanR")) +
  scale_fill_manual(values = c("#be9c00ff", "#00be70ff", "#00c1abff", "#00bbdaff", "#ff65acff")) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 12)) +
  labs(x = "Number of unique virulence factors \n co-occurring on ARG contigs per gOTU", 
       y = NULL)


f3e <- 
  ggplot(gotu.both[gotu.both$gotu.count > 5,], 
         aes(x = n.virulence.argclass, y = prevalence, fill = ARG.class, shape = family.sign)) +
  geom_vline(xintercept = 100, color = "black", size = 0.5, linetype="dotted") +
  geom_hline(yintercept = 0.75, color = "black", size = 0.5, linetype="dotted") +
  geom_point(size = 2.5,alpha = 0.7) +
  labs(x = "Number of unique virulence factors \n co-occurring on ARG contigs per gOTU", 
       y = "ARG class prevalence per gOTU") +
  scale_shape_manual(values=c(22, 21, 23)) +
  scale_fill_manual(values = c("#f8766dff", "#be9c00ff", "#00be70ff", "#00c1abff", "#00bbdaff", "grey", "#f962ddff", "#ff65acff")) +
  theme_bw() + 
  theme(legend.position ="none", text = element_text(size = 12))

gotu.both1 <- gotu.both[(gotu.both$n.virulence.argclass >= 100) & (gotu.both$prevalence >= 0.75) & (gotu.both$gotu.count >= 6),]

library("gridExtra")
library("cowplot")

pdf(file = "figures\\Figure_3_te.pdf", width = 11.7, height = 8.3)
ggdraw() +
  draw_plot(f3a, 0, .6, 0.33, .4) +
  draw_plot(f3b, 0.34, .6, 0.66, .4) +
  draw_plot(f3c, 0, 0, .32, .6) +
  draw_plot(f3d, .33, 0, .32, .6) +
  draw_plot(f3e, .67, 0, .32, .6) +
  draw_plot_label(c("A", "B", "C", "D", "E"), c(0, 0.33, 0, 0.33, 0.66), c(1, 1, 0.6, 0.6, 0.6), size = 16)
dev.off()

# ggdraw() +
#   draw_plot(f3a, 0, .67, 0.33, .33) +
#   draw_plot(f3b, 0.34, .67, 0.33, .33) +
#   draw_plot(f3c, 0.67, 0, .325, 1) +
#   draw_plot(f3d, 0, 0, .325, .65) +
#   draw_plot(f3e, .34, 0, .325, .65) +
#   draw_plot_label(c("A", "B", "C", "D", "E"), c(0, 0.33, 0, 0.33, 0.66), c(1, 1, 0.5, 0.5, 0.5), size = 16)

# get taxa of red dots

gotu.both1 <- gotu.both[gotu.both$highlight == 1,]

pdf(file = "figures\\SupplementaryFigure7.pdf", height = 11.7, width = 8.3)
fs7
dev.off()

