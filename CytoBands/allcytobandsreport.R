load('cytobands.DM.Rda')
cytoba<-as(cytobands,'data.frame')

expected.w.statistic<-(sum(normals)*sum(tumors))/2

wilcoxon.res<-apply(cytobands.methylation,1,function(row){
			meth.values<-jitter(row)
			w<-wilcox.test(meth.values[normals],meth.values[tumors])
			c(w$p.value,(w[['statistic']]<expected.w.statistic))
		})
cb.wilcoxon.p.values<-wilcoxon.res[1,]
cb.normals.are.less.methylated<-as.logical(wilcoxon.res[2,])

cb.wilcoxon.p.values.bonferroni<-p.adjust(cb.wilcoxon.p.values,'bonferroni')
cb.wilcoxon.p.values.fdr<-p.adjust(cb.wilcoxon.p.values,'fdr')
cytoba<-cbind(cytoba,
	'wilcoxon.p.values'=cb.wilcoxon.p.values,
	'is.hyper'=cb.normals.are.less.methylated,
	'bonferroni'=cb.wilcoxon.p.values.bonferroni,
	'fdr'=cb.wilcoxon.p.values.fdr
)



methmeans<-apply(cytobands.methylation,1,function(row){
			c(mean(row[which(normals)]),mean(row[which(tumors)]))
		})

norm.meth.rate<-methmeans[1,]/width(cytobands)
tumor.meth.rate<-methmeans[2,]/width(cytobands)

cytoba<-cbind(cytoba,'norm.rate'=norm.meth.rate,'tumor.rate'=tumor.meth.rate)


write.table(cytoba,file='cytobands.review.tsv',sep='\t',row.names=FALSE,col.names=TRUE)

