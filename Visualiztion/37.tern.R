# https://mp.weixin.qq.com/s?__biz=MzkwNTY1ODMwMw==&mid=2247487956&idx=2&sn=b77ee182876c9d2a24edd77265327d9c&chksm=c0f53b04f782b2120aae75a8e2e357c32b0e6d19262b4ea16500aa5090f77032670e52d15874&scene=21#wechat_redirect

rm(list=ls())

# load pkgs---------------------------------------------------------------------
# install.packages('ggtern')
library(ggplot2)
library(ggtern)
library(reshape2)
library(dplyr)
# load data---------------------------------------------------------------------
dat <- data.frame(
  a = c(20,20,50,80,100),
  b = c(20,30,50,10,0),
  stringsAsFactors = FALSE
) %>% 
  mutate(c = 100 - a - b)

dat <- dat / rowSums(dat)

dat <- dat %>% 
  mutate(group = letters[1:5] %>% factor,
         weight = c(1,1.2,0.6,2,3))
# Visualization-----------------------------------------------------------------
Palette <- c("#9999CC","#CC6666","#CC79A7","#D55E00", 
             "#0072B2","#F0E442","#009E73","#56B4E9",
             "#E69F00","#B2182B")

ggtern(dat, aes(x = a, y = b, z = c)) + 
  geom_mask() + 
  geom_point(aes(size = weight, color = group), alpha = 0.8) +
  scale_size(range = c(1, 10)) +
  scale_color_manual(values  = Palette) +
  guides(size="none") +
  guides(color = guide_legend(nrow =10,
                              override.aes = list(size = 6),
                              title.position = "top")) + # 调整图例位置
  theme_bw()+
  theme(legend.text = element_text(size=10,face = "bold"),
        legend.position = c(0.85,0.72))

