# https://en.wikipedia.org/wiki/Decision_curve_analysis
# Threshold probability is defined as the minimum probability of an event at which a decision-maker would take a given action, for instance, the probability of cancer at which a doctor would order a biopsy. 

library(dcurves)

head(df_binary)

gmodels::CrossTable(df_binary$cancer, df_binary$cancerpredmarker > 0)
gmodels::CrossTable(df_binary$cancer, df_binary$cancerpredmarker > .14)

# treat none
y = 0
# treat all
# TP / N - FP / N * (pt / (1 - pt))
(105 / 750) - (645 / 750) * (0  / (1 - 0)) # 0.14
(105 / 750) - (645 / 750) * (0.14 / (1 - 0.14)) # 0
# (105 / 750) - (645 / 750) * (0.5 / (1 - 0.5))

res <- dca(cancer ~ famhistory, df_binary)

dca(cancer ~ famhistory, df_binary) %>%
  plot(smooth = TRUE)




  
  
  
  


