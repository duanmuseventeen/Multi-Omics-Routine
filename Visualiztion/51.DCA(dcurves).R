require(dcurves)
require(ggplot2)
require(dplyr)

dat      <- readxl::read_excel("51.data for dca.xlsx")

dat4roc <- dat %>% 
  mutate(group = case_when(label == "AA" ~ 1, TRUE ~ 0)) %>% 
  filter(complete.cases(`SA`))
fit.sa <- glm(group ~ `SA`, data = dat4roc, family = binomial(link="logit"))

dat4roc$pre <- predict(fit.sa, type = "response")

dca_result <- dca(group ~ pre, 
                  data = dat4roc,
                  thresholds = seq(0, 0.99, by = 0.01),
                  label = list(pre = "SA"))

p_dca.sa <- dca_result %>%
  plot(smooth = TRUE) + # smooth = TRUE 使曲线更平滑
  ggplot2::labs(
    x = "Threshold Probability", 
    y = "Net Benefit"
  ) +
  scale_color_manual(values = c("gray90","gray90","#E64B35"), 
                     labels = c("Treat All", "Treat None", "Combination")) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = element_text(size = 14),
    legend.position = c(0.8, 0.8),
    panel.grid = element_blank()
  )
ggsave(p_dca.sa, filename = "DCA(SA).pdf", width=7, height=7, units="in")
