setwd("C:/D/R project/Multi-Omics-Routine/Biostatistic/")

load("mymetabolomics.Rdata")
# lm----------------------------------------------------------------------------
# method 1
dat.wide %>% 
  mutate(x = `Stearoyl-L carnitine`, y = `Oleoyl-L carnitine`) %>% 
  ggplot(aes(x = x, y = y))  +
  geom_point() + 
  geom_smooth(method="lm", 
              formula = y ~ rcs(x, 3),
              se = TRUE, colour="blue") 

# method 2
dd <- datadist(dat.wide)
options(datadist = "dd")
fit.lm <- lm(`Oleoyl-L carnitine` ~ rcs(`Stearoyl-L carnitine`), data = dat.wide)
# number of knots in rcs defaulting to 5

pre <- data.frame(
  `Stearoyl-L carnitine` = seq(min(dat.wide$`Stearoyl-L carnitine`), 
                               max(dat.wide$`Stearoyl-L carnitine`),
                               length.out = 500),
  check.names = FALSE
  )

# 使用模型进行预测，得到预测值和置信区间
pre$yhat<- predict(fit.lm, newdata = pre)
pre$se  <- predict(fit.lm, newdata = pre, se.fit = TRUE)$se.fit
pre$lci <- pre$yhat - 1.96 * pre$se
pre$uci <- pre$yhat + 1.96 * pre$se

# 使用 ggplot2 绘制 RCS 曲线
ggplot(dat.wide, aes(x = `Stearoyl-L carnitine`, y = `Oleoyl-L carnitine`)) +
  geom_point(alpha = 0.5) +
  geom_line(data = pre, aes(y = yhat), color = "red3", size = 1) +
  geom_ribbon(data = pre, aes(y = yhat, ymin = lci, ymax = uci), 
              fill = "red3", alpha = 0.2) +
  theme_minimal()
# glm---------------------------------------------------------------------------
dat4lr <- dat.wide %>%
  mutate(outcome = PCA, C18 = `Stearoyl-L carnitine`) %>% 
  dplyr::select(PCA, C18)

dd <- datadist(dat4lr)
options(datadist = "dd")

fit.lr <- lrm(PCA ~ rcs(C18, 3), data = dat4lr)

pre <- rms::Predict(fit.lr, C18, conf.int = 0.95, ref.zero = TRUE, fun = exp)

ggplot(pre, aes(x = C18)) +
  geom_line(aes(y = yhat), color = "red3", size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "red3", alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray80") +
  theme_minimal()
# cox---------------------------------------------------------------------------
