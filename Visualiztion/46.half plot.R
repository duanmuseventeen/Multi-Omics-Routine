# Ref
# https://ggplot2.tidyverse.org/reference/geom_point.html
# https://ggplot2.tidyverse.org/reference/layer_stats.html
# https://ggplot2.tidyverse.org/reference/geom_function.html

require(ggplot2)
require(dplyr)
require(ggsci)
require(ggpubr)
require(scales)
# devtools::install_github("JanCoUnchained/ggunchained")
require(ggunchained)
require(gghalves)

# ggunchained-------------------------------------------------------------------
tibble::tibble(
  symbol = rep(letters[1:8], each = 100),
  group = rep(c("A","B"),each = 50) %>% rep(8),
  value = rnorm(800)) %>% 
  ggplot(aes(x = symbol, y = value, fill = group))+
  geom_split_violin(trim = T, colour = NA)+
  geom_point(stat = 'summary',fun=mean,
             position = position_dodge(width = 0.9)) +
  stat_summary(fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = 'errorbar',color='black',
               width = 0.01,size = 0.5,
               position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#197EC099", "#FED43999")) +
  stat_compare_means(label = 'p.signif') +
  labs(xlab = "", ylab = "Relative expression") +
  theme_bw() +
  theme(text = element_text(size = 12))


tibble::tibble(
  symbol = rep(letters[1:8], each = 100),
  group = rep(c("A","B"),each = 50) %>% rep(8),
  value = rnorm(800)) %>% 
  ggplot(aes(x = symbol, y = value, fill = group))+
  geom_split_violin(trim = T, colour = NA)+
  stat_summary(geom = "point", fun = mean,
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = 'errorbar',color='black',
               width = 0.01,size = 0.5,
               position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#197EC099", "#FED43999")) +
  stat_compare_means(label = 'p.signif') +
  labs(xlab = "", ylab = "Relative expression") +
  theme_bw() +
  theme(text = element_text(size = 12))
# gghalves----------------------------------------------------------------------
dat <- tibble::tibble(
  symbol = rep(letters[1:8], each = 100),
  group = rep(c("A","B"),each = 50) %>% rep(8),
  value = rnorm(800)) 

ggplot()+
  geom_half_violin(
    data = dat %>% filter(group == "A"),
    aes(x = symbol, y = value),colour=NA,fill="#197EC099",side = "l"
  ) +
  geom_half_violin(
    data = dat %>% filter(group == "B"),
    aes(x = symbol, y = value),colour=NA,fill="#FED43999",side = "r"
  ) +
  geom_point(data = dat, aes(x = symbol, y = value, fill = group),
             stat = 'summary', fun = mean,
             position = position_dodge(width = 0.7)) +
  stat_summary(data = dat, aes(x = symbol, y = value, fill = group),
               fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = 'errorbar', # color='black',
               width = 0.1, size = 0.5,
               position = position_dodge(width = 0.7))+
  stat_compare_means(data = dat, aes(x = symbol, y = value, fill = group), 
                     label = 'p.signif') +
  labs(xlab = "", ylab = "Relative expression") +
  theme_bw() +
  theme(text = element_text(size = 12))
