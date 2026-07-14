mycalibration <- function(data, adjust = FALSE, B = 1000, title){
  require(rms)
  
  # function ----
  mycalibration.plot <- function(cal, title){
    df = data.frame(
      x                = cal[,"predy"],
      apparent         = cal[,"calibrated.orig"],
      `bias-corrected` = cal[,"calibrated.corrected"],
      check.names = FALSE
    ) 
    
    df = df %>% pivot_longer(cols = c("apparent", "bias-corrected")) %>% 
      mutate(name = factor(name, levels = c("bias-corrected", "apparent")))
    
    ggplot(df, aes(x = x, y = value, col = name, linetype = name)) +
      geom_abline(slope = 1, intercept = 0, col = "gray90", linetype = 2) +
      geom_line() +
      scale_color_manual(values = c( "#E64B35", "#4DBBD5")) +
      labs(x = "Predicted probability", 
           y = "Actual Probability", 
           title = paste0(title, "\nBrier Score: ", brier)) +
      coord_equal(xlim = c(0,1), ylim = c(0,1)) +
      theme_bw() +
      theme(panel.grid = element_blank(), text = element_text(size = 12),
            plot.title = element_text(hjust = 0.5),
            plot.margin = ggplot2::margin(10,10,10,10))
  }
  
  # calibration----
  if(!adjust){
    rmsfit = lrm(group_int ~ x, data = data, x = TRUE, y = TRUE)
  }else if(adjust){
    rmsfit = lrm(group_int ~ x + age + sex, data = data, x = TRUE, y = TRUE)    
  }
  
  cal = calibrate(rmsfit, method = "boot", B = B)
  
  val <- rms::validate(rmsfit, method = "boot", B = B)
  brier = val["B", "index.corrected"] %>% sprintf("%.2f", .)
  
  p = mycalibration.plot(cal = cal, title = title)
  
  # return ----
  res = list(
    fit  = rmsfit,
    cal  = cal,
    plot = p
    )
  return(res)
}

require(dplyr)
require(ggplot2)

cal <- list()
cal[[1]] <- mycalibration(jph, adjust = F, title = "A calibration: univariable")
cal[[2]] <- mycalibration(jph, adjust = T, title = "A calibration: multivariable")
cal[[3]] <- mycalibration(szh, adjust = F, title = "B calibration: univariable")
cal[[4]] <- mycalibration(szh, adjust = T, title = "B calibration: multivariable")

p_out <- 
  cal[[1]]$plot + cal[[2]]$plot  + 
  cal[[3]]$plot + cal[[4]]$plot + 
  plot_layout(ncol = 2, nrow = 2)

ggsave(plot = p_out, filename = paste0("Calibration",".pdf"), device = "pdf",width = 9, height = 9, dpi = 300)
