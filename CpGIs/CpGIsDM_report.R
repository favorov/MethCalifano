source('CpCIsDM.R')
columns<-c('id','space','start','end')
karyotype.stat<-cbind(CpGIs.with.methylation[,columns],'wilcoxon.p.value'=wilcoxon.p.values,'is.hyper'=normals.are.less.methylated)
print('No correction')
print(karyotype.stat[DM.CpGIslands,])
print('Bonferroni')
print(karyotype.stat[DM.CpGIslands.Bonferroni,])

