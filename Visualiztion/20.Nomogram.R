# Nomogram & Calibration for diagnosis
# Not run

require(regplot)
require(rms)

geom_calibrate <- function(
    cal, 
    colors = c("gray80","#0C56A0","#b9292b"),
    B,
    MAE,
    n
    ){
  
  data <- data.frame(
    x = cal[,'predy'],
    Apparent = cal[,'predy'],
    Ideal = cal[,'calibrated.orig'],
    `Biascorrected` = cal[,'calibrated.corrected'],
    stringsAsFactors = FALSE
  ) %>%
    reshape2::melt("x") %>%
    mutate(variable = factor(variable, levels = c("Apparent","Biascorrected","Ideal")))
  
  ggplot() +
    geom_line(data = data, aes(x = x, y = value, color = variable, linetype = variable)) +
    geom_rug(data = data.frame(predicted = attr(cal, "predicted")),
             aes(x = predicted), sides = 't', alpha = 1/2) +
    scale_x_continuous(
      name = paste0("Predicted Probability\n","B = ",B,
                    " repetitions, boot    Mean absolute error = ",MAE,
                    "    n =", n), 
      limits = c(0,1)) +
    scale_y_continuous(name = "Actual Probability", limits = c(0,1)) + 
    scale_color_manual(values = colors) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(text = element_text(size = 12), legend.title = element_blank(),
          plot.margin = ggplot2::margin(2,30,30,30),
          panel.grid = element_blank(), legend.background = element_blank(),
          legend.position = c(0.2, 0.8)
    )
}

dd = rms::datadist(dat);options(datadist="dd")

fit.nomo <- rms::lrm(group ~ x1 + x2 + x3, data = dat, maxit=1000, x = TRUE, y = TRUE)

nano <- rms::nomogram(
  fit.nomo, 
  fun = plogis, 
  lp = FALSE, 
  fun.at = c(0.05,seq(0.1,0.9,by = 0.1),0.95),
  funlabel = "Risk")

# Nomogram----------------------------------------------------------------------
plot(nano)

regplot(fit.nomo, 
        observation= NULL, 
        clickable=FALSE, 
        title = "",
        plots = c("spikes","spikes"), 
        spkcol = "#d9904c",
        points=TRUE, rank = NULL) -> myregplot

write.csv(myregplot[[6]], "D:/myregplot.csv")
# Calibration-------------------------------------------------------------------
# require(PredictABEL)
calib <- calibrate(fit.nomo, method = 'boot', B = 1000)

plot(calib, xlim = c(0,1.0), ylim = c(0,1.0))
# Divergence or singularity in 945 samples
# 
# n=98   Mean absolute error=0.013   Mean squared error=0.00019
# 0.9 Quantile of absolute error=0.018

geom_calibrate(calib, B = 1000, MAE = 0.013, n = 945)

