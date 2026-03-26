library(forestploter)

p_forest <- forest(
  data = df_for_plot[, c("study", "Hazard Ratio (95% CI)", " ", "P Value")],
  est = df_for_plot$est,
  lower = df_for_plot$lower,
  upper = df_for_plot$upper, 
  ci_column = 3,             
  ref_line = 1,              # HR = 1 的参考线
  arrow_lab = c("Low Risk", "High Risk"), # 底部箭头标注
  xlim = c(0, ceiling(max(df_for_plot$upper))), 
  ticks_at = c(0, 1, 2, 4)   # 设置刻度
)

print(p_forest)
