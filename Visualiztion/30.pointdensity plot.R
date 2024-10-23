# ATAC-Seq analysis reveals a widespread decrease of chromatin accessibility in age-related macular degeneration
# PMID: 29636475

require(dplyr)
require(stringr)
require(ggpointdensity)
require(ggplot2)
require(viridis)

Retina <- data.table::fread("GSE99287_Retina_ATACSeq_peak_counts.txt.gz")

Retina$AMD4_Retina_PerR_Rep1 <- log2(Retina$AMD4_Retina_PerR_Rep1 + 1)
Retina$AMD4_Retina_PerR_Rep2 <- log2(Retina$AMD4_Retina_PerR_Rep2 + 1)

ggplot(Retina, aes(x = AMD4_Retina_PerR_Rep1, y = AMD4_Retina_PerR_Rep2)) +
  annotate("text", x = 2.5, y = 13, label = "R = 0.98", hjust = 0, size = 6) +
  geom_pointdensity() +
  geom_abline(intercept = 0, slope = 1, color = "gray80", linetype = 2) +
  scale_x_continuous(limits = c(0,15), expand = c(0,0), name = "Replicate 1 from peripheral retina of AMD4") +
  scale_y_continuous(limits = c(0,15), expand = c(0,0), name = "Replicate 2 from peripheral\n retina of AMD4") +
  scale_color_viridis() +
  theme_classic() +
  theme(text = element_text(size = 20), panel.grid = element_blank(),
        plot.margin = ggplot2::margin(30,30,30,30),
        legend.title = element_blank(),legend.text = element_blank(),
        legend.position = c(0.8,0.4), legend.margin = ggplot2::margin(0,0,0,0),
        legend.ticks = element_blank(),legend.background = element_rect(linetype = 1, colour = "#555555"),
        axis.title.x = element_text(vjust = -4), axis.title.y = element_text(vjust = 4)
  )
