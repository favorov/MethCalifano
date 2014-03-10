source('CpGIsDM.R')
columns<-c('id','space','start','end')
CpGIs.stat<-cbind(CpGIs.with.methylation[,columns],'wilcoxon.p.value'=wilcoxon.p.values,'is.hyper'=normals.are.less.methylated)
print('Differential methylation (0.05) with Bonferroni')
print(CpGIs.stat[DM.CpGIslands.Bonferroni,])

