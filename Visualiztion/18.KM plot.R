# KM-Plot

# single run--------------------------------------------------------------------
dat1 <- readxl::read_excel("18.38395893-Fig4f.xlsx", sheet = "Sheet1")

mysurv <- survival::Surv(dat1$`OS (months)`, dat1$`OS_status (1=death)` == 1)

tmp <- dat1 %>% dplyr::select("OS (months)","OS_status (1=death)","Riskscore") %>%
  mutate(var_split = ifelse(Riskscore <= 2.1, 0, 1))

mysurv <- survival::Surv(tmp$`OS (months)`, tmp$`OS_status (1=death)` == 1)
data.sm <- data.frame(
  time = tmp$`OS (months)`,
  status = tmp$`OS_status (1=death)`,
  var_split = tmp$var_split
)
fit <- survfit(mysurv ~ var_split, data = data.sm)
survdiff(mysurv ~ var_split, data = tmp)$pvalue
ggsurvplot(fit, data = data.sm,
           palette = c("#53944f","#6ea5f4"),
           conf.int = F, conf.int.style = "step",
           risk.table = FALSE,
           legend = c(0.15, 0.1),
           ylim = c(0.3, 1),
           pval = TRUE, pval.method = F, pval.coord = c(0,0.5),
           legend.labs = c("Low risk","High risk"),
           tables.height = 0.25,
           pval.size = 6, pval.method.size = 5,
           font.x = c(20),  font.y = c(20),
           ggtheme =
             theme_bw() +
             theme(text = element_text(size = 20),
                   legend.title = element_blank(),
                   plot.margin = ggplot2::margin(30,30,30,30),
                   panel.grid = element_blank(),
                   legend.background = element_blank(),
                   axis.title.x = element_text(vjust = -4),
                   axis.title.y = element_text(vjust = 4))) +
  xlab("Time (months)") +
  ylab("OS probability")

# batch-------------------------------------------------------------------------
# Not run
# calculate p-value of log-rank-------------------------------------------------
report <- data.frame(
  `Gene Symbol` = Genes,
  bin = rep(NA, n),
  tri = rep(NA, n),
  qua = rep(NA, n),
  row.names = Genes,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
for (i in report$`Gene Symbol`) {
  for(j in c(2, 3, 4)){
    if(i %in% colnames(dat)){
      tmp <- KMplot(kmdata = dat, var = i, split = j)
      mysurv <- survival::Surv(tmp$bcr_follow_up_m, tmp$bcr == 1)
      data.sm <- data.frame(
        time = tmp$bcr_follow_up_m,
        status = tmp$bcr,
        var_split = tmp$var_split
      )
      
      if(j == 2){
        report$bin[report$`Gene Symbol` == i] <- survdiff(mysurv ~ var_split, data = tmp)$pvalue
      }else if(j == 3){
        report$tri[report$`Gene Symbol` == i] <- survdiff(mysurv ~ var_split, data = tmp)$pvalue
      }else if(j == 4){
        report$qua[report$`Gene Symbol` == i] <- survdiff(mysurv ~ var_split, data = tmp)$pvalue
      }
    }
  }
}

# calculate p-value of Cox regression-------------------------------------------
report <- data.frame(
  `GENE SYMBOL` = Genes,
  B = NA,
  se = NA,
  p = NA,
  HR = NA,
  lci = NA,
  uci = NA,
  row.names = Genes,
  check.names = F,
  stringsAsFactors = F
)
mysurv <- Surv(dat$bcr_follow_up_m, dat$bcr == 1)
for(i in 1:nrow(report)){
  if(report$`GENE SYMBOL`[i] %in% colnames(dat)){
    fit.cph.uni <- coxph(
      as.formula(paste0("mysurv ~ `", report$`GENE SYMBOL`[i],"`")), 
      data = dat)
    report[i,2:7] <- c(
      summary(fit.cph.uni)$coefficients[1,1], #B
      summary(fit.cph.uni)$coefficients[1,3], # se 
      summary(fit.cph.uni)$coefficients[1,5], # p
      summary(fit.cph.uni)$coefficients[1,2],# HR 
      summary(fit.cph.uni)$conf.int[1,3], # lci
      summary(fit.cph.uni)$conf.int[1,4] # uci
    )
  }
}

report <- report %>% 
  mutate(`95% CI` = paste0(round(HR,2)," (",round(lci,2),"-",round(uci,2),")"))
# batch for KM-plot-------------------------------------------------------------
KMplot <- function(kmdata, var, split = 2){
  kmdata <- kmdata %>% dplyr::select("bcr_follow_up_m","bcr",var)
  colnames(kmdata)[3] <- "var"
  if(split == 2){
    kmdata <- kmdata %>% 
      mutate(var_split = ifelse(var <= median(var), 0, 1))
  }else if(split == 3){
    kmdata <- kmdata %>% 
      mutate(var_split = ifelse(var <= quantile(var,.33), 0,
                                ifelse(var <= quantile(var,.67),1,2))) %>%
      filter(var_split %in% c(0,2))
  }else if(split == 4){
    kmdata <- kmdata %>% 
      mutate(var_split = ifelse(var <= quantile(var,.25), 0,
                                ifelse(var <= quantile(var,.5),1,
                                       ifelse(var <= quantile(var,.75),2,3)))) %>%
      filter(var_split %in% c(0,3))
  }
  
  return(kmdata)
}

pList <- list();n <- 1;splitn <- 2

for (i in targetGenes) {
  if(i %in% colnames(dat)){
    tmp <- KMplot(kmdata = dat, var = i, split = splitn)
    mysurv <- survival::Surv(tmp$bcr_follow_up_m, tmp$bcr == 1)
    data.sm <- data.frame(
      time = tmp$bcr_follow_up_m,
      status = tmp$bcr,
      var_split = tmp$var_split
    )
    pList[[n]] <- ggsurvplot(fit, data = data.sm,
                             palette = c("#d9904c","#394d54"),
                             conf.int = TRUE, conf.int.style = "step",
                             risk.table = FALSE,
                             legend = c(0.75, 0.25),
                             legend.title = i,
                             legend.labs = paste0(c("Low ","High "), i),
                             pval = TRUE, pval.method = TRUE,
                             tables.height = 0.25,
                             pval.size = 5, pval.method.size = 5,
                             font.x = c(5),  font.y = c(5)
    )
    
    # surv_diff <- survdiff(mysurv ~ var_split, data = tmp)
    # print(surv_diff)
    
    print(n)
    n <- n + 1
    }
}
p_out <- arrange_ggsurvplots(pList, print = FALSE, ncol = 8, nrow = 8)
ggsave(filename = "KM-plots.pdf", plot = p_out, device = "pdf", path = "D:/", 
       width = 16, height = 16, dpi = 300, units = "in")
