# correlation matrix plot

require(dplyr)
require(GGally)

dat1 <- readxl::read_excel("27.correlation matrix plot.xlsx", sheet = "1") 
dat2 <- readxl::read_excel("27.correlation matrix plot.xlsx", sheet = "2") 
dat3 <- readxl::read_excel("27.correlation matrix plot.xlsx", sheet = "3") 
dat4 <- readxl::read_excel("27.correlation matrix plot.xlsx", sheet = "4") 
dat5 <- readxl::read_excel("27.correlation matrix plot.xlsx", sheet = "5") 

ggpairs(dat1, columns = 2:8, title = "Technical replicates: CLGS")
ggpairs(dat2, columns = 2:8, title = "Technical replicates: IFRV")
ggpairs(dat3, columns = 2:8, title = "Technical replicates: VGTE")
ggpairs(dat4, columns = 2:8, title = "Technical replicates: MVRN")
ggpairs(dat5, columns = 2:8, title = "Technical replicates: NREC")

# OPTIMAL-----------------------------------------------------------------------
# https://zhuanlan.zhihu.com/p/506277796

# Not Run
p <- ggpairs(dat.ally, columns = 2:8, title = "",
        mapping = ggplot2::aes(color = group),
        upper = list(continuous = wrap("cor",method = "spearman")), # 更改相关性系数，?ggally_cor，查看函数参数，默认pearson。
        lower = list(continuous = "smooth", # ggally_smooth()
                     #combo = "facethist", 
                     #discrete = "facetbar", 
                     na ="na"),
        # lower = list(continuous = wrap("smooth",se = TRUE)), # 默认lm
        # lower = list(continuous = wrap("smooth",method = "loess",se = TRUE)), # 默认lm
        diag = list(continuous = "densityDiag", # 'densityDiag', 'barDiag', 'blankDiag'可选
                    #discrete = "barDiag", # 'barDiag', 'blankDiag'可选 
                    #diag =NULL # 不显示diag部分
                    na = "naDiag")
        ) + 
  scale_color_manual(values = c("#f67152","#fcb460","#aadfe5","#6e77b8"))

ggsave(plot = p, filename = "correlation matrix.pdf", device = "pdf",width = 10, height = 10, dpi = 300)
