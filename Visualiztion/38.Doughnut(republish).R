# https://mp.weixin.qq.com/s/_9h4Yf234zHaOPzh8h4OTg

rm(list=ls())

# load pkgs---------------------------------------------------------------------
library(ggplot2)
library(ggforce)
library(RColorBrewer)
# load data---------------------------------------------------------------------
dat <- data.frame(
  ratio = c(76.3, 13.8, 5.2, 3.5, 1.2),
  kingdom = c("Bacteria", "Eukaryota", "Viruses", "Unassigned", "Archaea")
)

brewer.pal.info
# maxcolors category colorblind
# BrBG            11      div       TRUE
# PiYG            11      div       TRUE
# PRGn            11      div       TRUE
# PuOr            11      div       TRUE
# RdBu            11      div       TRUE
# RdGy            11      div      FALSE
# RdYlBu          11      div       TRUE
# RdYlGn          11      div      FALSE
# Spectral        11      div      FALSE
# Accent           8     qual      FALSE
# Dark2            8     qual       TRUE
# Paired          12     qual       TRUE
# Pastel1          9     qual      FALSE
# Pastel2          8     qual      FALSE
# Set1             9     qual      FALSE
# Set2             8     qual       TRUE
# Set3            12     qual      FALSE
# Blues            9      seq       TRUE
# BuGn             9      seq       TRUE
# BuPu             9      seq       TRUE
# GnBu             9      seq       TRUE
# Greens           9      seq       TRUE
# Greys            9      seq       TRUE
# Oranges          9      seq       TRUE
# OrRd             9      seq       TRUE
# PuBu             9      seq       TRUE
# PuBuGn           9      seq       TRUE
# PuRd             9      seq       TRUE
# Purples          9      seq       TRUE
# RdPu             9      seq       TRUE
# Reds             9      seq       TRUE
# YlGn             9      seq       TRUE
# YlGnBu           9      seq       TRUE
# YlOrBr           9      seq       TRUE
# YlOrRd           9      seq       TRUE

display.brewer.all(type="qual") 

# Visualization-----------------------------------------------------------------
# pie
p1 <- ggplot(dat, aes(x = "", y = ratio, fill = kingdom)) +
  geom_bar(width = 1, stat = "identity", color = "white") + # 设置颜色为白色
  coord_polar("y", start = 0) + # 设置极坐标系
  geom_text(aes(label = sprintf("%1.1f%%", (ratio/sum(ratio))*100)),
            position = position_stack(vjust = 0.5), # 垂直居中文本
            color = "black", size = 5) + # 格式化标签为百分比
  scale_fill_manual(values = brewer.pal(5, "Set2")) + # 手动设置颜色
  theme_void() # 去除图形背景和网格线

p1

# doughnut
p2 <- ggplot() +  
  geom_arc_bar(data = dat,  
               stat = "pie",  
               aes(x0 = 0, y0 = 0, r0 = 1, r = 2, # r0 = 0就是饼图 
                   amount = ratio, fill = kingdom)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  labs(x = "", y = "") +
  scale_fill_manual(values = brewer.pal(5, "Set2")) +  
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.ticks = element_blank(),   
        axis.text.y = element_blank(),  
        axis.text.x = element_blank(),  
        legend.title = element_blank(),   
        panel.border = element_blank(),  
        panel.background = element_blank())

p2









