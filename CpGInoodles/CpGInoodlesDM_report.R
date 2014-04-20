source('CpGIsDM.R')

cat('Differential methylation (0.05) with Bonferroni\n')
print.data.frame(DM.CpGIs.stat,row.names=FALSE)


