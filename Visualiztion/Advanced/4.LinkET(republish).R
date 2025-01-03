## LinkET-----------------------------------------------------------------------
require(linkET)
require(tibble)

## load data--------------------------------------------------------------------
dat <- tibble(
  id = letters,
  a  = runif(26),
  b1 = runif(26),
  b2 = runif(26),  
  b3 = runif(26),
)
## calculate cor----------------------------------------------------------------
mantel <- data.frame(
  main = rep("a", each = (ncol(dat) - 2)),
  variate  = colnames(dat)[-c(1:2)],
  stringsAsFactors = F
) %>% 
  mutate(r = NA, p = NA)
for(i in unique(mantel$main)){
  for (j in unique(mantel$variate)) {
    print(c(i,j))
    tmp <- cor.test(dat %>% dplyr::select(i) %>% unlist,
                    dat %>% dplyr::select(j) %>% unlist)
    mantel$r[mantel$main == i & mantel$variate == j] <- tmp$estimate
    mantel$p[mantel$main == i & mantel$variate == j] <- tmp$p.value
  }
}
mantel <- mantel %>%
  mutate(r = abs(r)) %>%
  mutate(rd = cut(r, breaks = c(0, 0.2, 0.4, Inf), 
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

## graphics---------------------------------------------------------------------
p <-qcorrplot(correlate(dat[,-c(1:2)]), 
              type = "upper", diag = FALSE, method = "spearman") +
  geom_square() +
  geom_couple(aes(alpha = r, colour = pd, size = rd), # key
              data = mantel,
              curvature = nice_curvature()) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")), 
                       limits = c(0,1)) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Spearman's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Spearman's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Spearman's r", order = 3),
         alpha = guide_legend(title = "Spearman's r", order = 3))

ggsave(plot = p, filename = "LinkET.pdf", device = "pdf",width = 10, height = 10, dpi = 300)

