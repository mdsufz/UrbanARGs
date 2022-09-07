# This script takes the CheckM data dn plots a figure of the completeness and
# contamination of all MAGs

# packages
library(data.table)
library(ggplot2)
library(patchwork)
library(viridis)
library(gridExtra)
library(dplyr)

# load checkm data
checkm <- fread("data\\checkm.tsv", sep = "\t")
gotu.final <- fread("data\\final_groups_qual.tsv", sep = "\t")
gotu.final.95 <- gotu.final[gotu.final$ani == 95,]

# create plot of mag quality
checkm <- as.data.frame(checkm)
checkm$quality <- "Low"
checkm$quality[(checkm$Completeness > 50) & (checkm$Contamination < 10)] <- "Medium"
checkm$quality[(checkm$Completeness > 90) & (checkm$Contamination < 5)] <- "High"
checkm$quality.score <-checkm$Completeness - 5*checkm$Contamination

# summarize HQ and MQ MAGs
mag.quality <- as.data.frame(checkm %>%
                               group_by(quality) %>%
                               summarize(count = n_distinct(mag)) %>% 
                               mutate(percent = round(prop.table(count) * 100, 1)))

mag.scatter <- ggplot(checkm, aes(x = Contamination, y = Completeness, color = quality.score)) +
  geom_point(size = 3, show.legend = FALSE) +
  scale_color_viridis(option = "viridis") +
  scale_y_continuous(limits = c(50, 100)) +
  scale_x_continuous(limits = c(0, 10)) +
  theme_bw() +
  theme(text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Contamination (%)") + ylab("Completeness (%)") +
  annotation_custom(tableGrob(mag.quality, rows = NULL, theme = ttheme_minimal()), 
                    xmin=5, xmax=10, ymin=50, ymax=70)

mag.completeness.hist <- ggplot(checkm, aes(x = Completeness, fill = cut(Completeness, 100))) + 
  geom_histogram(bins = 100, show.legend = FALSE) + 
  scale_fill_viridis(option = "viridis", discrete = T) +
  theme_void() +
  scale_x_continuous(limits = c(50, 100)) +
  coord_flip()

mag.contamination.hist <- ggplot(checkm, aes(x = Contamination, fill = cut(Contamination, 100))) + 
  geom_histogram(bins = 200, show.legend = FALSE) + 
  scale_fill_viridis(option = "viridis", discrete = T, direction = -1) +
  scale_x_continuous(limits = c(0, 10)) +
  theme_void()

pdf(file = paste0("figures\\completeness_contamination_mag.pdf"),
    height = 8.27, width = 11.69)

mag.contamination.hist + plot_spacer() + mag.scatter + mag.completeness.hist + 
  plot_layout(
    ncol = 2, 
    nrow = 2, 
    widths = c(4, 1),
    heights = c(1, 4)
  ) 

dev.off()


# for gOTU species only
gotu.checkm <- checkm[which(checkm$mag %in% gotu.final.95$bin),]

# summarize HQ and MQ MAGs
gotu.quality <- as.data.frame(gotu.checkm %>%
  group_by(quality) %>%
  summarize(count = n_distinct(mag)) %>% 
  mutate(percent = round(prop.table(count) * 100, 1)))

final.clusters.scatter <- ggplot(gotu.checkm, aes(x = Contamination, y = Completeness, color = quality.score)) +
  geom_point(size = 3, show.legend = FALSE) +
  scale_color_viridis(option = "viridis") +
  scale_y_continuous(limits = c(50, 100)) +
  scale_x_continuous(limits = c(0, 10)) +
  theme_bw() +
  theme(text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Contamination (%)") + ylab("Completeness (%)") +
  annotation_custom(tableGrob(gotu.quality, rows = NULL, theme = ttheme_minimal()), 
                    xmin=5, xmax=10, ymin=50, ymax=70)

final.clusters.completeness.hist <- ggplot(gotu.checkm, aes(x = Completeness, fill = cut(Completeness, 100))) + 
  geom_histogram(bins = 100, show.legend = FALSE) + 
  scale_fill_viridis(option = "viridis", discrete = T) +
  theme_void() +
  scale_x_continuous(limits = c(50, 100)) +
  coord_flip()

final.clusters.contamination.hist <- ggplot(gotu.checkm, aes(x = Contamination, fill = cut(Contamination, 100))) + 
  geom_histogram(bins = 200, show.legend = FALSE) + 
  scale_fill_viridis(option = "viridis", discrete = T, direction = -1) +
  scale_x_continuous(limits = c(0, 10)) +
  theme_void()

pdf(file = paste0("figures\\completeness_contamination_gotu.pdf"),
    height = 8.27, width = 11.69)

final.clusters.contamination.hist + plot_spacer() + final.clusters.scatter + final.clusters.completeness.hist + 
  plot_layout(
    ncol = 2, 
    nrow = 2, 
    widths = c(4, 1),
    heights = c(1, 4)
  ) 

dev.off()
