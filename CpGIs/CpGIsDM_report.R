source('CpGIsDM.R')

print('Differential methylation (0.05) with Bonferroni')
print.data.frame(DM.CpGIs.stat,row.names=FALSE)


