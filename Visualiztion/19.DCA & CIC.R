# Decision Curve Analysis (DCA) & Clinical Impact Curve (CIC)
# Not run

require(rmda)

set.seed(1011)
model1 <- decision_curve(group ~ x1, data = dat, thresholds=seq(0,1,by=0.05), family = binomial(link='logit'), study.design ='case-control', bootstraps = 1000)
model2 <- decision_curve(group ~ x1 + x2, data = dat, thresholds=seq(0,1,by=0.05), family = binomial(link='logit'), study.design ='case-control', bootstraps = 1000)
model3 <- decision_curve(group ~ x1 + x2 + x3, data = dat, thresholds=seq(0,1,by=0.05), family = binomial(link='logit'), study.design ='case-control', bootstraps = 1000)

plot_decision_curve(list(model1, model2, model3),
                    curve.names = c("Model1", "Model2","Model3"),
                    col=c("#4DAF4A","#3A84C1", "#E41A1C"), 
                    lty=c(1,1,1), 
                    confidence.intervals = FALSE, 
                    cost.benefit.axis = FALSE, 
                    standardize = FALSE, 
                    legend.position = "bottomleft")


rmda::plot_clinical_impact(model1, population.size = 1000,
                           cost.benefit.axis = T,
                           n.cost.benefits = 8,
                           col = c('red','blue'),
                           confidence.intervals = T,
                           ylim = c(0,1000),
                           legend.position ='topright')
rmda::plot_clinical_impact(model2, population.size = 1000,
                           cost.benefit.axis = T,
                           n.cost.benefits = 8,
                           col = c('red','blue'),
                           confidence.intervals = T,
                           ylim = c(0,1000),
                           legend.position ='topright')
rmda::plot_clinical_impact(model3, population.size = 1000,
                           cost.benefit.axis = T,
                           n.cost.benefits = 8,
                           col = c('red','blue'),
                           confidence.intervals = T,
                           ylim = c(0,1000),
                           legend.position ='topright')

