# Ref:
# https://mp.weixin.qq.com/s/KA6UW4q-BzmnIdpqvoI5lQ

rm(list = ls())
# Source Data-------------------------------------------------------------------
# Post-translational modifications orchestrate the intrinsic signaling bias of GPR52
# pmid: 40087539
# Figure 1A left panel

# Note: IC50 can be calculated using the same function of EC50

# load pkgs---------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(magrittr)
library(cowplot)
library(ggprism)
library(scales)
# load data---------------------------------------------------------------------
dat <- readxl::read_excel("43.40087539 Figure1A.xlsx")
dat <- dat %>% 
  mutate(dose = case_when(
    `log [GLP1 (7–36)] M` == "-∞" ~ 0,
    `log [GLP1 (7–36)] M` == "-12" ~ exp(-12),
    `log [GLP1 (7–36)] M` == "-11" ~ exp(-11),
    `log [GLP1 (7–36)] M` == "-10" ~ exp(-10),
    `log [GLP1 (7–36)] M` == "-9" ~ exp(-9),
    `log [GLP1 (7–36)] M` == "-8" ~ exp(-8),
    `log [GLP1 (7–36)] M` == "-7" ~ exp(-7),
    `log [GLP1 (7–36)] M` == "-6" ~ exp(-6)),
    log_dose = case_when(
      `log [GLP1 (7–36)] M` == "-∞" ~ -Inf,
      `log [GLP1 (7–36)] M` == "-12" ~ -12,
      `log [GLP1 (7–36)] M` == "-11" ~ -11,
      `log [GLP1 (7–36)] M` == "-10" ~ -10,
      `log [GLP1 (7–36)] M` == "-9" ~ -9,
      `log [GLP1 (7–36)] M` == "-8" ~ -8,
      `log [GLP1 (7–36)] M` == "-7" ~ -7,
      `log [GLP1 (7–36)] M` == "-6" ~ -6))

# geom_smooth-------------------------------------------------------------------
cols <- RColorBrewer::brewer.pal(4,'Set2') #from RColorBrewer

ggplot(dat, aes(x = dose, y = `cAMP response (% agonist max)`, 
                color = group, shape = group)) +
  geom_point(size = 4) +
  geom_smooth(method = "nls", formula = y ~ SSfpl(x, A, B, xmid, scal), se = FALSE) +
  scale_x_log10(guide = "prism_minor",
                breaks = trans_breaks(log10, function(x) 10^x),
                labels = trans_format(log10, math_format(10^.x)),
                minor_breaks = rep(1:9, 4)*(10^rep(-2:-5, each = 9)),
                expand = c(0, 0))+
  scale_color_manual(values = c("#66C2A5", "#FC8D62")) + 
  guides(color = "none", shape = "none") +
  theme_bw() +
  theme(text = element_text(size = 16),
        panel.grid = element_blank())


ggplot(dat, aes(x = log(dose), y = `cAMP response (% agonist max)`, 
                color = group, shape = group)) +
  geom_point(size = 4) +
  geom_smooth(method = "nls", formula = y ~ SSfpl(x, A, B, C, D), se = FALSE) +
  scale_color_manual(values = c("#66C2A5", "#FC8D62")) + 
  guides(color = "none", shape = "none") +
  theme_bw() +
  theme(text = element_text(size = 16),
        panel.grid = element_blank())










