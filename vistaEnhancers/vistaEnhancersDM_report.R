source('vistaEnhancersDM.R')
cat('Differential methylation (0.05) with Bonferroni\n')
print.data.frame(DM.vistaEnhancers.stat,row.names=FALSE)

