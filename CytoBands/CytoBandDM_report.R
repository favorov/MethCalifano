source('CytoBandDM.R')
columns<-c('id','space','start','end')
karyotype.stat<-cbind(karyotype.with.methylation[,columns],'p.value'=wilcoxon.p.values,'bonferroni'=wilcoxon.p.values.bonferroni,'fdr'=wilcoxon.p.values.fdr,'is.hyper'=normals.are.less.methylated)
print('No correction')
karyotype.stat.DM<-karyotype.stat[DM.cytobands,]
karyotype.stat.DM<-karyotype.stat.DM[order(as.integer(substring(karyotype.stat.DM$space,4)),karyotype.stat.DM$start),]
print.data.frame(karyotype.stat.DM[,c(columns,'p.value','is.hyper')],row.names=FALSE)
print('Benjamini-Hoochberg')
print(karyotype.stat[DM.cytobands.fdr,c(columns,'p.value','fdr','is.hyper')])
print('Bonferroni')
print(karyotype.stat[DM.cytobands.bonferroni,c(columns,'p.value','bonferroni','is.hyper')])

