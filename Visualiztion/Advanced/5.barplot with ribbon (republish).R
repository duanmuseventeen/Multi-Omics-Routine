# optimization of stack bar plot

# Reference
# https://zhuanlan.zhihu.com/p/551933995

# load pkgs---------------------------------------------------------------------
library(ggplot2)
library(ggalluvial)

# load data---------------------------------------------------------------------
dat <- data.frame(
  group = sample(c("control","model","treatment"), 100, replace = TRUE),
  variate = rep(c("a","b"), each = 50))
# bar plot----------------------------------------------------------------------
dat %>% 
  group_by(group, variate) %>% 
  mutate(n = n()) %>% 
  ggplot(aes(x = variate, fill = group))+
  geom_bar(width = 0.5)
# bar plot with ribbon----------------------------------------------------------
dat %>% 
  mutate(value = 1) %>% 
  group_by(group, variate) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = variate, y = n, fill = group, stratum = group, alluvium = group))+
  geom_stratum(width = 0.5, color='white') +
  geom_alluvium(alpha = 0.5,
                width = 0.5,
                curve_type = "linear") +
  labs(x = "", y = "Percent") +
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position="none",
        panel.grid=element_blank())
