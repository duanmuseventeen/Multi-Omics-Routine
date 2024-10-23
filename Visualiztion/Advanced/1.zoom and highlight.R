# zoom and highlight

require(dplyr)
require(ggplot2)
require(ggforce)

dat <- readxl::read_excel("7.38347143-Source Data Fig2.xlsx", sheet = "2e")

dat <- dat %>% 
  filter(CNA == "Amplification") %>%
  mutate(CBCGA = CBCGA * 100,
         TCGA = TCGA * 100,
         col = case_when(Gene %in% c("MYC","MCL1","ERBB2","BTG2") ~ 1,
                         TRUE ~ 0),
         lab = case_when(Gene %in% c("MYC","MCL1","ERBB2","BTG2") ~ Gene,
                         TRUE ~ NA),
         col = factor(col)) 
# Repeat figure in publication--------------------------------------------------
p1 <- dat %>% 
  ggplot(aes(x = CBCGA, y = TCGA, colour = col)) +
  annotate(geom = "rect", 
           xmin = 0, xmax = 25,
           ymin = 0, ymax = 25,
           linetype = 2, 
           linewidth = 2, 
           col = "gray20", fill = NA) +
  annotate(geom = "segment", 
           x = c(25,25,100), xend = c(100,100,100),
           y = c(0,25,25), yend = c(0,25,100),
           linetype = 2, 
           linewidth = 2, 
           col = "gray20") +
  geom_point(size = 4) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 2, col = "gray50") +
  scale_color_manual(values = c("gray80", "red3")) +
  scale_x_continuous(limits = c(0,100)) +
  scale_y_continuous(limits = c(0,100)) +
  guides(col = "none") +
  labs(x = "Prevalence in CBCGA (%)", y = "Amplification in\nTCGA white (%)") +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(text = element_text(size = 16),
        panel.grid = element_blank())

p2 <- dat %>% 
  ggplot(aes(x = CBCGA, y = TCGA, colour = col, label = lab)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(col = "black", size = 6) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 2, col = "gray50") +
  scale_color_manual(values = c("gray80", "red3")) +
  scale_x_continuous(limits = c(0,25)) +
  scale_y_continuous(limits = c(0,25)) +
  guides(col = "none") +
  labs(x = "Prevalence in CBCGA (%)", y = "Amplification in\nTCGA white (%)") +
  theme_bw() +
  theme(text = element_text(size = 16),
        panel.grid = element_blank())

p1 + p2

# Revision::zoom by ggforce-----------------------------------------------------
dat %>% 
  ggplot(aes(x = CBCGA, y = TCGA, colour = col)) +
  geom_point(size = 4) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 2, col = "gray50") +
  scale_color_manual(values = c("gray80", "red3")) +
  scale_x_continuous(limits = c(0,100)) +
  scale_y_continuous(limits = c(0,100)) +
  guides(col = "none") +
  labs(x = "Prevalence in CBCGA (%)", y = "Amplification in\nTCGA white (%)") +
  ggforce::facet_zoom(xlim = c(0,25), ylim = c(0,25), zoom.size = 1, show.area = TRUE) +
  theme_bw() +
  theme(text = element_text(size = 16),
        panel.grid = element_blank())
# Revision::zoom by annotate----------------------------------------------------
p <- ggplot(dat, aes(x = CBCGA, y = TCGA, colour = col, label = lab)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(col = "black", size = 6) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 2, col = "gray50") +
  scale_color_manual(values = c("gray80", "red3")) +
  scale_x_continuous(limits = c(0,25)) +
  scale_y_continuous(limits = c(0,25)) +
  guides(col = "none") +
  labs(x = "", y = "") +
  theme_bw() +
  theme(text = element_text(size = 16),
        panel.grid = element_blank())

dat %>% 
  ggplot(aes(x = CBCGA, y = TCGA, colour = col)) +
  geom_point(size = 4) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 2, col = "gray50") +
  annotation_custom(grob = ggplot2::ggplotGrob(p), 
                    xmin = 25, xmax = 75, ymin = 25, ymax = 75) +
  scale_color_manual(values = c("gray80", "red3")) +
  scale_x_continuous(limits = c(0,100)) +
  scale_y_continuous(limits = c(0,100)) +
  guides(col = "none") +
  labs(x = "Prevalence in CBCGA (%)", y = "Amplification in\nTCGA white (%)") +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(text = element_text(size = 16),
        panel.grid = element_blank())


