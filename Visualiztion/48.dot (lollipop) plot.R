require(stringr)
require(dplyr)
require(ggplot2)

dat <- readxl::read_excel("data.xlsx") %>% 
  arrange(log2FC)
dat <- dat %>% 
  mutate(Name = factor(Name, levels = Name))
p <- dat %>% 
  ggplot(aes(x = log2FC, y = Name, size = -log10(`p-value`), colour = log2FC)) +
  geom_vline(xintercept = 0, col = "gray90", linetype = 2) + 
  geom_segment(data = dat, aes(xend = 0, yend = Name), col = "grey80", linewidth = 0.5) +
  geom_point() +
  scale_color_gradient2(low = "#0C56A0", mid = "#ffffff", high = "#b9292b", midpoint = 0) +
  labs(x = "log2FC(COPD vs VC)", y = "") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(p, filename = "dotplot.pdf", width=5, height=8, units="in")
