# Integrated multiomic profiling of breast cancer in the Chinese population reveals patient stratification and therapeutic vulnerabilities
# PMID: 38347143

require(ggtext)
require(ggplot2)

Fig2c2 <-
  readxl::read_excel("7.38347143-Source Data Fig2.xlsx", sheet = "2c") %>% 
  mutate(Gene = factor(Gene, levels = c(
    'SF3B1',
    'PTEN',
    'GATA3',
    'MAP3K1',
    'FOXA1',
    'KMT2C',
    'PIK3CA',
    'TP53',
    'CBFB',
    'CDKN1B',
    'AKT1'))) %>% 
  ggplot(aes(x = Ability, y = Gene, col = Pathway)) +
  geom_segment(aes(x = Ability - s.e., y = Gene, xend = Ability + s.e., yend = Gene)) +
  geom_point(size = 3) +
  guides(color = guide_legend(nrow = 2)) +
  scale_x_reverse() +
  labs(y = "", x = "", color = "") +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.position = "bottom",
        axis.text.y = element_markdown(face = "italic"),
        axis.line.y = element_blank(), axis.ticks.y = element_blank())