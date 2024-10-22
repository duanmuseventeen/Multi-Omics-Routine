# ROC curve

require(pROC)
require(ggplot2)

dat <- readxl::read_excel("17.38395893-Fig4b.xlsx")

fit <- glm(Status ~ `Predicted risk score`, data = dat, family = binomial(link="logit"))

pre <- predict(fit, type='response')

modelroc <- roc(dat$Status, pre)

plot(modelroc)

tmp <- data.frame(
  Specificity = modelroc$specificities,
  Sensitivity = modelroc$sensitivities,
  stringsAsFactors = F)

ggplot(tmp, aes(Specificity, Sensitivity, col = "Lasso: AUROC=0.832\n95%Cl: 0.697 - 0.951")) +
  geom_abline(slope = 1, intercept = 1, col = "#8d8d93", linetype = 2, linewidth = 2) +
  geom_path(linewidth = 2) +
  scale_x_reverse(expand = c(0,0), n.breaks = 6) +
  scale_y_continuous(expand = c(0,0), n.breaks = 6) +
  scale_color_manual(values = "#22962d") +
  labs(x = "False positive rate", y = "True positive rate") +
  coord_fixed(ratio = 1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.margin = ggplot2::margin(30,30,30,30),
        legend.title = element_blank(),
        legend.background = element_rect(color = "gray90"),
        text = element_text(size = 16),
        legend.position = c(0.6,0.1))

modelroc
# Call:
#   roc.default(response = dat$Status, predictor = pre)
# 
# Data: pre in 50 controls (dat$Status 0) < 10 cases (dat$Status 1).
# Area under the curve: 0.832
# 
# ci(modelroc)
# 95% CI: 0.6802-0.9838 (DeLong)