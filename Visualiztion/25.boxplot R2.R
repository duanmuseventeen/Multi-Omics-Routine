# boxplot

require(ggplot2)
require(ggpubr)
require(dplyr)

dat <- readxl::read_excel("25.38347143-Source Data Extended Fig6.xlsx", sheet = "ED_6b")

dat %>% 
  mutate(x = factor(`Intrinsic subtype`, 
                    levels = c("LumA","LumB","HER2","Basal","Normal")),
         group = factor(`TME phenotype`, 
                        levels = c("Cold","Moderate","Hot"))) %>% 
  ggplot(aes(x, `Immune Signature`, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(#jitter.height=0.75, 
                                              jitter.width = 0.2, 
                                              dodge.width = 0.8)) +
  annotate(geom = "segment", 
           x = seq(0.7,4.7,1), xend = seq(1.3,5.3,1),
           y = rep(1.65,5), yend = rep(1.65,5)) +
  stat_compare_means(aes(group = `TME phenotype`), 
                     method = "kruskal.test", 
                     label = "p.signif", 
                     hide.ns = F, size = 6) +
  scale_color_manual(values = c("#6184a5","#dbad34","#aa1722")) +
  labs(x = "", y = "Immune Signature", color = "TME phenotype") +
  theme_classic() +
  theme(text = element_text(size = 16, color = "black"), 
        legend.position = "top", legend.title = element_text(face = "bold"))
     

