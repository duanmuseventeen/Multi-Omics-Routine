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
             position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c("#197EC099", "#FED43999")) +
  stat_summary(fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = 'errorbar',color='black',
               width = 0.01,size = 0.5,
               position = position_dodge(width = 0.9)) +
  labs(xlab = "", ylab = "Relative expression") +
  theme_bw() +
  theme(text = element_text(size = 12))

# ggunchained-------------------------------------------------------------------
tibble::tibble(
  symbol = rep(letters[1:8], each = 100),
  group = rep(c("A","B"),each = 50) %>% rep(8),
  value = rnorm(800)) %>% 
  ggplot(aes(x = symbol, y = value, fill = group))+
  geom_split_violin(trim = T, colour = NA)+
  geom_point(stat = 'summary',fun=mean,
             position = position_dodge(width = 0.9))+
  scale_fill_manual(values = c("#197EC099", "#FED43999")) +
  stat_summary(fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = 'errorbar',color='black',
               width = 0.01,size = 0.5,
               position = position_dodge(width = 0.9)) +
  labs(xlab = "", ylab = "Relative expression") +
  theme_bw() +
  theme(text = element_text(size = 12))










