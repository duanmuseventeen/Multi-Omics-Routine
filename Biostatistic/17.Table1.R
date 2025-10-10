rm(list = ls())

require(stringr)
require(dplyr)

checksig <- function(res){
  print(paste0("t.pval < 0.05: ",length(res$ions[res$t.pval < 0.05])))
  print(paste0("FDRq < 0.05: ",length(res$ions[res$FDRq < 0.05])))
  print(paste0("FDRq < 0.05 and abs(log2FC) > log2(1.5): ",length(res$ions[res$FDRq < 0.05 & abs(log2(res$FC)) > log2(1.5)])))
  print(paste0("FDRq < 0.05 and abs(log2FC) > 1: ",length(res$ions[res$FDRq < 0.05 & abs(log2(res$FC)) > 1])))
}

myout.IQR <- function(x){
  y <- x
  x[y < (quantile(y, .25) - 1.5 * IQR(y))] <- (quantile(y, .25) - 1.5 * IQR(y))
  x[y > (quantile(y, .75) + 1.5 * IQR(y))] <- (quantile(y, .75) + 1.5 * IQR(y))
  return(x)
}

geom_volcano <- function(dat, pos.num = 1.5, neg.num = -1.5, 
                         xmin = -2, xmax = 2, title, pshape = 21){
  require(ggplot2)
  require(ggrepel)
  
  dat <- dat %>% 
    mutate(log2FoldChange = log2(FC),
           padj = FDRq,
           SYMBOL = ions) %>% 
    mutate(color = case_when(log2FoldChange > 0 & padj < 0.05 ~ 2,
                             log2FoldChange < 0 & padj < 0.05 ~ 1,
                             TRUE ~ 0)) %>% 
    mutate(label = factor(color, levels = c(0,1,2), labels = c("Non-Sig","Down","Up")),
           color = factor(color))
  dat_up <- dat %>% 
    filter(label == "Up") %>% 
    arrange(desc(log2FoldChange))
  dat_up <- dat_up[1:ifelse(nrow(dat_up) >= 10, 10, nrow(dat_up)),]
  if(is.na(dat_up$SYMBOL[1]) & nrow(dat_up) == 1){
    dat_up[1,] <- c(NA,0,0,0,0,1,1,NA,NA)
    dat_up$log2FoldChange <- dat_up$log2FoldChange %>% as.numeric
    dat_up$padj <- dat_up$padj %>% as.numeric
  }
  dat_down <- dat %>% 
    filter(label == "Down") %>% 
    arrange(log2FoldChange)
  dat_down <- dat_down[1:ifelse(nrow(dat_down) >= 10, 10, nrow(dat_down)),]
  if(is.na(dat_down$SYMBOL[1]) & nrow(dat_down) == 1){
    dat_down[1,] <- c("",0,0,0,0,1,1,NA,NA)
    dat_down$log2FoldChange <- dat_down$log2FoldChange %>% as.numeric
    dat_down$padj <- dat_down$padj %>% as.numeric
  }
  
  ggplot(dat, aes(x = log2FoldChange, y = -log10(padj), fill = label, label = SYMBOL)) +
    geom_point(shape = pshape) +
    geom_vline(xintercept = c(0), color = "gray80", linetype = 2) +
    geom_hline(yintercept = c(1.30103), color = "gray80", linetype = 2) +
    ylab(expression(-log[10]~(adj.~P~value))) +
    xlab("Log2(Fold Change)") +
    labs(fill = "") +
    scale_x_continuous(limits = c(xmin, xmax)) +
    scale_fill_manual(values = c("gray90", "blue3", "red3")) +
    scale_size_continuous(range = c(0.1, 4)) +
    geom_text_repel(
      data = dat_up,
      color = "red3",
      size = 5,
      nudge_x = pos.num - as.numeric(dat_up$log2FoldChange),
      segment.size=0.3,
      segment.color="grey",
      direction="y",
      hjust= 0,
      max.overlaps = Inf) +
    geom_text_repel(
      data= dat_down,
      color="blue3",
      size=5,
      nudge_x = neg.num - as.numeric(dat_down$log2FoldChange),
      segment.size = 0.3,
      segment.color = "grey",
      direction="y",
      hjust= 1,
      max.overlaps = Inf) +
    labs(title = title) + 
    theme_minimal() +
    theme(legend.position = "right",
          legend.title.align = 0, # left align
          legend.title = element_text(margin = margin(t = 15, unit = "pt")) # add more space on top of legend titles
    ) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(size=20),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18),
          aspect.ratio = 1/1.2, panel.grid.major = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    theme(plot.title = element_text(hjust = 0.5, size=20))
}
# load data---------------------------------------------------------------------
dat <- readxl::read_excel("jacc.rawdat.xlsx")
# del metabolite without identification-----------------------------------------
colnames(dat)

ions <- str_extract(colnames(dat), "^M[0-9]{2,4}T[0-9]{1,3}$")
ions <- ions[!is.na(ions)]
ions

length(ions)
# [1] 1844
dat.metab <- dat %>% dplyr::select(-all_of(ions))

dim(dat)
# [1] 2020 2033
dim(dat.metab)
# [1] 2020  189
# filter samples----------------------------------------------------------------
meta <- readxl::read_excel("jacc.meta.xlsx") %>% 
  mutate(group = str_remove_all(ID, "-[0-9]+$"))
meta.new <- meta %>% 
  filter(group %in% c("NCA", "NOCA", "AMI"))

dim(meta.new)
# [1] 1298   25

meta.new <- meta.new %>% 
  mutate(group2 = case_when(
    Cre <= 80 & group == "NCA" ~ 0, # 正常组
    Cre > 80 & group %in% c("NCA", "NOCA") ~ 1, # 肾损伤组
    Cre > 80 & group == "AMI" ~ 2, # 心肾共损组
    TRUE ~ 3
  )) %>% 
  filter(group2 != 3)

dim(meta)
# [1] 2020   24
dim(meta.new)
# [1] 671  26

table(meta.new$group2)
# 0   1   2 
# 166 250 255 
# 基线信息----------------------------------------------------------------------
table(meta.new$group2)
# 0   1   2 
# 166 250 255

table1contious <- function(dat = meta.new, myvar){
  var <- dat %>% dplyr::select(myvar) %>% unlist
  nc <- var[dat$group2 == 0]
  ckd<- var[dat$group2 == 1]
  ckdami <- var[dat$group2 == 2]
  
  fit <- lm(var ~ dat$group2)
  cat(paste0(
    myvar, ":",
    "\nNC: ", mean(nc, na.rm = TRUE) %>% format(digits = 5), "±", sd(nc, na.rm = TRUE) %>% format(digits = 5),
    "\nCKD: ", mean(ckd, na.rm = TRUE) %>% format(digits = 5), "±", sd(ckd, na.rm = TRUE) %>% format(digits = 5),
    "\nCKD + AMI: ", mean(ckdami, na.rm = TRUE) %>% format(digits = 5), "±", sd(ckdami, na.rm = TRUE) %>% format(digits = 5),
    "\nP for trend: ", summary(fit)$coefficients[2,4] %>% format(digits = 3, scientific = TRUE)
  ))
}
table1contious(myvar = "age")
table1contious(myvar = "lvef")
table1contious(myvar = "sp")
table1contious(myvar = "dp")
table1contious(myvar = "hr")
table1contious(myvar = "HbA1c")
table1contious(myvar = "glucose")
table1contious(myvar = "TG")
table1contious(myvar = "TC")
table1contious(myvar = "HDL-C")
table1contious(myvar = "LDL-C")
table1contious(myvar = "BUN")
table1contious(myvar = "Cre")

table1cate <- function(dat = meta.new, myvar, pos){
  var <- dat %>% dplyr::select(myvar) %>% unlist
  dat <- dat %>% 
    mutate(group2 = factor(group2, labels = c("NC","CKD","CKD&AMI")))
  res <- table(var, dat$group2) %>% as.data.frame %>% 
    group_by(Var2) %>% 
    mutate(Sum = sum(Freq),
           Prop = format(Freq / sum(Freq) * 100, digits = 4)) %>% 
    mutate(Prop = paste0(Freq," (",Prop,"%)")) %>% 
    filter(var == pos)
  
  trend <- prop.trend.test(
    x = res$Freq, 
    n = res$Sum, 
    score = c(0, 1, 2) # 赋值
  )
  trend$p.value
  
  return(list(`result` = res,`p-value for trend` = trend$p.value))
  
  # https://zhuanlan.zhihu.com/p/450066515
  # tmp <- table(var, dat$group2) %>% as.data.frame %>% 
  #   tidyr::pivot_wider(names_from = "var", values_from = "Freq")
  # trend <- matrix(c(tmp$`1`, tmp$`2`), byrow = FALSE, ncol = 2) %>% 
  #   DescTools::CochranArmitageTest()
  # trend$p.value
  # 两种结果一致
}
table1cate(myvar = "gender", pos = 2) # 性别（1为男性、2为女性）
table1cate(myvar = "hp", pos = 1) # 高血压史（1为有，2为无）
table1cate(myvar = "dm", pos = 1) # DM：糖尿病史（1为有，2为无）
table1cate(myvar = "smoking", pos = 1) # 吸烟史（1为有，2为无）

table1cate(myvar = "drinking", pos = 1) # 饮酒史（1为有，2为无）
