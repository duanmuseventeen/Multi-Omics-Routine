# data is a df including only QC samples with rowname, in which row means sample 
# and column means feature

qc.cdf <- function(dat){
  dat <- dat %>% 
    t %>% cor
  dat[row(dat) <= col(dat)] <- NA
  dat %>% 
    as.data.frame %>% 
    tibble::rownames_to_column("A") %>% 
    as.data.frame %>% 
    tidyr::pivot_longer(cols = !A, names_to = "B", values_to = "cor") %>% 
    filter(complete.cases(.)) %>% 
  ggplot(aes(cor)) +
  stat_ecdf() +
  scale_x_continuous(limits = c(-1,1)) +
  labs(x = "Correlation", y = "CDF") +
  theme_bw() +
  theme(text = element_text(size = 20),
        plot.margin =  ggplot2::margin(30,30,30,30),
        legend.background = element_rect(linetype = 1, color = "black"),
        legend.position = "right",
        axis.text = element_text(size= 16, color = "black"))
  } 
