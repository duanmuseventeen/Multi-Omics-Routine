# Ref:
# Theory
# https://www.graphpad.com/guides/prism/latest/curve-fitting/reg_the_ec50.htm

# https://github.com/DoseResponse/drc/blob/master/

# https://rstudio-pubs-static.s3.amazonaws.com/788707_37e78b45ea234e778c901227fba475ed.html
# https://www.jingege.wang/2023/08/09/%e5%89%82%e9%87%8f%e7%9b%b8%e5%ba%94%e5%85%b3%e7%b3%bb%e6%9b%b2%e7%ba%bf%e5%8f%8aic50%e8%ae%a1%e7%ae%97%ef%bc%9agraphpad-r%e8%af%ad%e8%a8%80/
# https://cran.r-project.org/web/packages/ec50estimator/vignettes/how_to_use.html

rm(list = ls())
# Source Data-------------------------------------------------------------------
# Post-translational modifications orchestrate the intrinsic signaling bias of GPR52
# pmid: 40087539
# Figure 1A left panel

# Note: IC50 can be calculated using the same function of EC50

# load pkgs---------------------------------------------------------------------
library(ec50estimator)
library(drc)
library(ggplot2)
library(dplyr)
library(magrittr)
library(cowplot)
library(ggridges)
# load data---------------------------------------------------------------------
dat <- readxl::read_excel("43.40087539 Figure1A.xlsx")
dat <- dat %>% 
  mutate(dose = case_when(
    `log [GLP1 (7–36)] M` == "-∞" ~ 0,
    `log [GLP1 (7–36)] M` == "-12" ~ exp(-12),
    `log [GLP1 (7–36)] M` == "-11" ~ exp(-11),
    `log [GLP1 (7–36)] M` == "-10" ~ exp(-10),
    `log [GLP1 (7–36)] M` == "-9" ~ exp(-9),
    `log [GLP1 (7–36)] M` == "-8" ~ exp(-8),
    `log [GLP1 (7–36)] M` == "-7" ~ exp(-7),
    `log [GLP1 (7–36)] M` == "-6" ~ exp(-6)),
    log_dose = case_when(
      `log [GLP1 (7–36)] M` == "-∞" ~ -Inf,
      `log [GLP1 (7–36)] M` == "-12" ~ -12,
      `log [GLP1 (7–36)] M` == "-11" ~ -11,
      `log [GLP1 (7–36)] M` == "-10" ~ -10,
      `log [GLP1 (7–36)] M` == "-9" ~ -9,
      `log [GLP1 (7–36)] M` == "-8" ~ -8,
      `log [GLP1 (7–36)] M` == "-7" ~ -7,
      `log [GLP1 (7–36)] M` == "-6" ~ -6))

# observe the data
dat %>% 
  ggplot(aes(log_dose, `cAMP response (% agonist max)`, color = group))+
  geom_jitter(width = 0.1)+
  theme_light()

# Calculate EC50 with ec50estimator---------------------------------------------
ec50 = estimate_EC50(`cAMP response (% agonist max)` ~ dose,
                     data = dat,
                     isolate_col = "group", 
                     strata_col = NULL,
                     interval = "delta",
                     fct = drc::LL.3())
ec50
# ID strata     Estimate   Std..Error         Lower        Upper
# 1 GLP-1R        4.417589e-05 2.539329e-06  3.889507e-05 4.945671e-05
# 2 vector        1.431688e-04 9.777365e-05 -6.016260e-05 3.465003e-04
# Calculate EC50 with drc-------------------------------------------------------
dat$group # Isolate
dat$dose # Concentration 

#When trying models out with non-transformed "conc", I kept getting fail to converge error messages.
#I need to log-transform to conc + 1 - this is how data is most often displayed in dose response curves.

#Log-logistic is the most commonly used model for dose response curves.
#However, my dataset is small (<15-20 data points per concentration level) so it may not 
#be appropriate here per Ritz et al. (2015) outlining this package.

#Trying initial four-parameter model
vector <- dat %>% filter(group == "vector")
ec50_vector <- drm(`cAMP response (% agonist max)` ~ dose, data = vector, fct = LL.4())
GLP_1R <- dat %>% filter(group == "GLP-1R")
ec50_GLP_1R <- drm(`cAMP response (% agonist max)` ~ dose, data = GLP_1R, fct = LL.4())

#### Selecting appropriate model------------------------------------------------
# Use mselect() to select the best-fitting model
#Here, I try some of the most commonly suggested in the example for mselect()
mselect(ec50_vector, list(LL.2(), LL.3(), LL.4(), W1.2(), W1.3(), baro5()))
# logLik       IC Lack of fit  Res var
# LL.4  -54.02379 118.0476  0.00377524 6.337138
# LL.2         NA       NA          NA       NA
# LL.3         NA       NA          NA       NA
# LL.4         NA       NA          NA       NA
# W1.2         NA       NA          NA       NA
# W1.3         NA       NA          NA       NA
# baro5        NA       NA          NA       NA
mselect(ec50_GLP_1R, list(LL.2(), LL.3(), LL.4(), W1.2(), W1.3(), baro5()))
# logLik       IC  Lack of fit  Res var
# LL.3  -72.13658 152.2732 2.350550e-06 27.30408
# LL.4  -72.13332 154.2666 1.070648e-06 28.66150
# LL.4  -72.13332 154.2666 1.070648e-06 28.66150
# LL.2         NA       NA           NA       NA
# W1.2         NA       NA           NA       NA
# W1.3         NA       NA           NA       NA
# baro5        NA       NA           NA       NA

#For Akaike's information criterion (AIC) and the residual standard error: 
#the smaller the better and for lack-of-fit test (against a one-way ANOVA model): 
#the larger (the p-value) the better. Note that the residual standard error 
#is only available for continuous dose-response data.
#Log likelihood values cannot be used for comparison unless the models are nested.

#Based off of AIC, W1.3 model is best which is a three-parameter Weibull model.
#(I also tried W1.4 and it wasn't better)

#Let's do it
ec50_GLP_1R <- drm(`cAMP response (% agonist max)` ~ dose, data = GLP_1R, fct = LL.3())

summary(ec50_GLP_1R)
# Parameter estimates:
#   
#   Estimate  Std. Error t-value   p-value    
# b:(Intercept) -2.4032e+00  3.6063e-01 -6.6641 1.350e-06 ***
#   d:(Intercept)  9.9259e+01  1.9280e+00 51.4834 < 2.2e-16 ***
#   e:(Intercept)  4.4176e-05  2.5393e-06 17.3967 5.975e-14 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error:
#   
#   5.225331 (21 degrees of freedom)

#Model fitted: Weibull (type 1) with lower limit at 0 (3 parms)
#note: Weibull is an asymmetric model type, while something like log-logistic is symmetric.

#The three parameters in this model are upper asymptote exp(c) 
#(not given here), lower asymptote- exp(d), slope - exp(b),
#and exp(e) - the inflection point, which is the relative EC50!

summary(ec50_GLP_1R)$coefficients["e:(Intercept)","Estimate"]
# [1] 4.417589e-05
# is the relative EC50. # # # # 

#We want both absolute and relative EC50 values.
#Absolute EC50 is where 50% of the maximum response occurs; 
#this is in line with how FRAC defines EC50
#and so it is what I will use as my phenotyping dose on media.
#As mentioned above, relative EC50 is the inflection point on a dose-response curve.

#From Noel et al. 2018: "The absolute EC50 was estimated by fitting percent relative 
#growth against the log concentration using a four-, three-, and two-parameter 
#log-logistic model but specifying type = “absolute” within the “ED” function of 
#“drc” (Ritz and Streibig 2005)."
#default in the ED() argument is relative

ED(ec50_GLP_1R, respLev = c(50), type = "absolute")
# Estimated effective doses
# 
# Estimate Std. Error
# e:1:50 4.4451e-05 2.5758e-06
# 4.4451e-05
# ppm is the absolute EC50. # # # # 

# CI
confint(ec50_GLP_1R, level = 0.95)
# 2.5 %        97.5 %
# b:(Intercept) -3.153204e+00 -1.653275e+00
# d:(Intercept)  9.524915e+01  1.032680e+02
# e:(Intercept)  3.889507e-05  4.945671e-05

#### Plotting the model with absolute and relative EC50 as labelled data points----
#Use my drm to predict y values for the relative EC50 (we know absolute y is 50)
x_value <- c(summary(ec50_GLP_1R)$coefficients["e:(Intercept)","Estimate"])
predicted_y <- predict(ec50_GLP_1R, newdata = data.frame(Dose = x_value))

# Make plot based on the drm model
plot(ec50_GLP_1R, 
     xlab = expression(bold("[GLP1 (7–36)] M")),
     ylab = expression(bold("cAMP response (% agonist max)")),
     lwd = 2, type = "n", pch = 19, xaxt = "n", yaxt = "n", legend = F, cex.axis = 1.2)
points(ED(ec50_GLP_1R, respLev = c(50), type = "absolute")[1,1], 50, pch = 19, cex = 1.5)
points(x_value, predicted_y, pch = 15, cex = 1.5)
text(ED(ec50_GLP_1R, respLev = c(50), type = "absolute")[1,1], 50, labels=expression(bold("Absolute EC50")), cex= 1, adj = c(-0.1, 1))
text(x_value, predicted_y, labels=expression(bold("Relative EC50")), cex= 1, adj = c(-0.1, -0.9))
title("Dose Response Curve")

cols <- RColorBrewer::brewer.pal(4,'Set2') #from RColorBrewer

plot(ec50_GLP_1R,
     type= "all",
     col=cols[1], #auto color selection for curves
     pch=16,
     lwd=1,
     cex.axis=0.8,
     # legend = TRUE,
     xlab= "[Compound] (nM)",
     ylab = "Intensity Sum")
plot(ec50_GLP_1R,
     col = cols[1],
     add=TRUE,
     type='confidence')

#### For multi group------------------------------------------------------------
ec50 <- drm(`cAMP response (% agonist max)` ~ dose, 
            curveid = group, data = dat, fct = LL.3())

2^summary(ec50)$coefficients[c("e:GLP-1R","e:vector"),"Estimate"]

ED(ec50, respLev = c(50), type = "absolute")

#### logDose--------------------------------------------------------------------
ec50 <- drm(`cAMP response (% agonist max)` ~ log_dose, 
            curveid = group, data = dat, logDose = exp(1), fct = LL.3())

2^summary(ec50)$coefficients[c("e:GLP-1R","e:vector"),"Estimate"]
# e:GLP-1R e:vector 
# 1.000031 1.000100 

ED(ec50, respLev = c(50), type = "absolute")
# Estimate Std. Error
# e:GLP-1R:50 4.4451e-05 2.0330e-06
# e:vector:50        NaN        NaN

x_value <- c(summary(ec50)$coefficients[c("e:GLP-1R","e:vector"),"Estimate"]) %>% log
predicted_y <- predict(ec50, newdata = 
                         data.frame(Log2Dose = x_value, 
                                    group = c("GLP-1R","vector")))
absolute_x  <- ED(ec50, respLev = c(50), type = "absolute")[1:2,1] %>% na.omit %>% log
plot(ec50, 
     xlab = expression(bold("Log [GLP1 (7–36)] M")),
     ylab = expression(bold("cAMP response (% agonist max)")),
     lwd = 2, type = "n", pch = 19, xaxt = "n", yaxt = "n", legend = F, cex.axis = 1.2)
points(absolute_x, rep(50, length(absolute_y)), pch = 19, cex = 1.5)
points(x_value, predicted_y, pch = 15, cex = 1.5)
text(absolute_x, rep(50, length(absolute_y)), labels=expression(bold("Absolute EC50")), cex= 1, adj = c(-0.1, 1))
text(x_value, predicted_y, labels=expression(bold("Relative EC50")), cex= 1, adj = c(-0.1, -0.9))
title("Dose Response Curve")

# 以上计算结果与Prism Graphpad 9计算结果相同
# 但原文结果无法重复
#### Visualization--------------------------------------------------------------
# method1: 
# method1.1:
ggplot() + 
  geom_jitter(data = dat, aes(x = dose, y = `cAMP response (% agonist max)`, 
                              colour = `cAMP response (% agonist max)`)) +
  geom_line(data = data.frame(
    dose = seq(range(dat$dose)[1],range(dat$dose)[2],length.out = 1000),
    stringsAsFactors = FALSE
  ) %>% 
    mutate(response = ec50_GLP_1R$curve[[1]](dose)),
  aes(x = dose, y = response, colour = response)) +
  scale_x_log10() + 
  scale_color_viridis_c(option = "B") +
  theme_classic() +
  theme(panel.grid = element_blank())

# method1.2:
line <- plot(ec50_GLP_1R,
             type= "all",
             pch=16,
             lwd=1,
             cex.axis=0.8,
             legend = TRUE,
             xlab= "[Compound] (nM)",
             ylab = "Intensity Sum")

se <- predict(ec50_GLP_1R, 
              newdata = data.frame(
                dose = seq(range(dat$dose)[1],range(dat$dose)[2],length.out = 1000)),
              interval = "confidence",
              level = 0.95)[, c("Lower", "Upper")] %>% 
  as.data.frame %>% 
  mutate(x = seq(range(dat$dose)[1],range(dat$dose)[2],length.out = 1000))

ggplot() + 
  geom_jitter(data = dat %>% filter(group == 'GLP-1R'),
              aes(x = log(dose), y = `cAMP response (% agonist max)`,
                              colour = `cAMP response (% agonist max)`)) +
  geom_line(data = line, aes(x = log(dose), y = `1`), colour = "steelblue") +
  geom_ribbon(data = se, aes(x = log(x), ymin = Lower, ymax = Upper), 
              fill = "steelblue", alpha = 0.2) +
  theme_classic() +
  theme(panel.grid = element_blank())
# method2: ---------------------------------------------------------------------
# https://rstudio-pubs-static.s3.amazonaws.com/788707_37e78b45ea234e778c901227fba475ed.html
library(RColorBrewer)

Finally, let’s plot the curves side- by side

cols <- brewer.pal(4,'Set2') #from RColorBrewer

par(mfrow=c(1,2)) #(n row, n col) 

for (i in 1:curve_cnt){ #note curves list index starts at 1
  DR = filter(DF_norm, Weeks==weeks[i,1]) 
  DR.mi <- drm(AggObj_IntSum ~ Concentration, 
               data= DR,
               robust = 'mean', #non-robust least squares estimation ("mean")
               fct = LL.4(names = c("Hill slope", "Min", "Max", "EC50")))
  plot(DR.mi,
       type= "all",
       col=cols[i], #auto color selection for curves
       pch=16,
       lwd=1,
       cex.axis=0.8,
       xlim = c(1e-4, 20),
       ylim = c(0, 5E+09),
       legend = TRUE,
       legendText = (paste(toString(weeks[i,1]), 'weeks')), 
       xlab= "[Compound] (nM)",
       ylab = "Intensity Sum")
  plot(DR.mi,
       col = cols[i],
       xlim = c(1e-4, 20),
       ylim = c(0, 5E+09),
       add=TRUE,
       type='confidence')}


As suspected, the curves does not reach a steady plateau. Therefore we can’t quantitatively compare best-fit curve EC50 and top values. However, there is a clear trend- the longer timecourse (7 weeks) has a stronger response in terms of both potency of compound (smaller EC50) and higher max (non-existent plateau). This can be easily visualized if we plot the suves on top of each other this time.

par(mfrow=c(1,1)) #(n row, n col) 

for (i in 1:curve_cnt){ #note curves list index starts at 1
  DR = filter(DF_norm, Weeks==weeks[i,1]) 
  DR.mi <- drm(AggObj_IntSum ~ Concentration, 
               data= DR,
               robust = 'mean', #non-robust least squares estimation ("mean")
               fct = LL.4(names = c("Hill slope", "Min", "Max", "EC50")))
  plot(DR.mi,
       type= "all",
       col=cols[i], #auto color selection for curves
       pch=16,
       lwd=1,
       cex.axis=0.8,
       xlim = c(1e-4, 20),
       ylim = c(0, 5E+09),
       legend = FALSE,
       xlab= "[Compound] (nM)",
       ylab = "Intensity Sum")
  plot(DR.mi,
       col = cols[i],
       xlim = c(1e-4, 20),
       ylim = c(0, 5E+09),
       add=TRUE,
       type='confidence')
  par(new=TRUE)} #Add new=TRUE arg to overlay graphs
