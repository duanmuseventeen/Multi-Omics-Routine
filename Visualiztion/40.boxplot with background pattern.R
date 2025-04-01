library(dplyr)
library(ggplot2)
library(stringr)

dat1 <- readxl::read_excel("Metagenome.xlsx")

cols <- colnames(dat1)[-1]
topn <- 10 # topn must be even

dat1$mean <- dat1[,-1] %>% rowMeans 
dat1 <- dat1 %>% 
  arrange(desc(mean)) %>% 
  mutate(order = row_number()) %>% 
  filter(order %in% c(1:topn))
lev <- dat1$ID

dat1 <- dat1 %>% 
  dplyr::select(-all_of(c('mean', 'order'))) %>% 
  tidyr::pivot_longer(cols = cols, names_to = "sample", values_to = "value") %>% 
  mutate(group = str_remove(sample, "[0-9]{1}$"),
         ID = factor(ID, levels = lev))

ggplot(dat1, aes(x = ID, y = value, color = group)) +
  annotate(geom = "rect", 
           xmin = seq(0.5, (topn - 2 + 0.5), 2), xmax = seq(1.5, (topn - 1 + 0.5), 2), 
           ymin = rep(0, length(seq(0.5, (topn - 2 + 0.5), 2))), 
           ymax = rep(25, length(seq(0.5, (topn - 2 + 0.5), 2))),
           fill = "skyblue", alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(#jitter.height=0.75, 
    jitter.width = 0.2, 
    dodge.width = 0.8)) +
  labs(x = "", y = "Abundance") +
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,25), expand = c(0,0)) + 
  theme_bw() +
  theme(text = element_text(size = 16, color = "black"), 
        panel.grid = element_blank(), 
        plot.margin = ggplot2::margin(30,30,30,30),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top", legend.title = element_text(face = "bold"))
