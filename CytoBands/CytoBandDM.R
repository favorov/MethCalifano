if (!suppressWarnings(require('differential.coverage')))
{
	if (!suppressWarnings(require('devtools')))
	{
		source("http://bioconductor.org/biocLite.R")
		biocLite("devtools")
		library("devtools")
	}
	install_github('favorov/differential.coverage')
	#load_all('../../../../differential.coverage/')
	library('differential.coverage')
}


cytobands.DM.loaded<-FALSE
# we can the whole thing to cytobands.DM.Rda
if(file.exists('cytobands.DM.Rda'))
	if ('cytobands.methylation' %in% load('cytobands.DM.Rda'))
		if (class(cytobands.methylation)=='dgCMatrix')
			cytobands.DM.loaded<-TRUE

if (!cytobands.DM.loaded)
{
	source('../common/read_clinical.R')
	#Clinical prepared.
	source('../common/prepare_beds_and_contrast.R')
	#bedfiles are ready
	cytobands<-get.cytoband.ranges()
	#cytobands ready
	SNPchip.version=package.version('SNPchip')
	cytobands.methylation<-CountCoverageOfNoodles(cytobands,bedfilenames,DNAids)
	#it is folder with bed files
	tests.number<-dim(cytobands.methylation)[1]

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
	cytobands.DM.statistics<-data.frame(
		'wilcoxon.p.values'=cb.wilcoxon.p.values,
		'is.hyper'=cb.normals.are.less.methylated,
		'bonferroni'=cb.wilcoxon.p.values.bonferroni,
		'fdr'=cb.wilcoxon.p.values.fdr
	)
	
	#DM.cytobands<-which(wilcoxon.p.values<=0.05)
	#DM.cytobands.bonferroni<-which(wilcoxon.p.values.bonferroni<=0.05)
	#DM.cytobands.fdr<-which(wilcoxon.p.values.fdr<=0.1)

	save(file='cytobands.DM.Rda',list=c('cytobands','cytobands.methylation','Clinical','clinFile','clinFileName','bedsinfolder','bed.used','tumors','normals','DNAids','SNPchip.version'))
}
