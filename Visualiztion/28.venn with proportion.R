# Venn plot, whose area are related to the number

require(eulerr)

dat1 <- readxl::read_excel("28.venn with proportion.xlsx", sheet = "1")
dat2 <- readxl::read_excel("28.venn with proportion.xlsx", sheet = "2")
dat3 <- readxl::read_excel("28.venn with proportion.xlsx", sheet = "3")

euler(
  list(`CNA-mRNA` = dat1$Symbol[dat1$`FDR\r\n(CNA-mRNA)` < 0.05],
       `CNA-Protein` = dat1$Symbol[dat1$`FDR\r\n(CNA-protein)` < 0.05]),
  shape = "circle") |>
  plot(
    quantities = list(type = c("counts", "percent"), cex=1),          
    labels = list(cex=1),                 
    edges  = list(col = "white", lex = 1), 
    fills  = list(fill = c("blue3","red3"), alpha=0.7),
    legend = list(side = "right")
  ) |>
  ggplotify::as.ggplot()
