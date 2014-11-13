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


karyotype.with.methylation.loaded<-FALSE
# we can the whole thing to karyotype.with.methylation.Rda
if(file.exists('karyotype.with.methylation.Rda'))
	if ('karyotype.with.methylation' %in% load('karyotype.with.methylation.Rda'))
		if (class(karyotype.with.methylation)=='data.frame')
			karyotype.with.methylation.loaded<-TRUE

if (!karyotype.with.methylation.loaded)
{
	source('../common/read_clinical.R')
	#Clinical prepared.
	source('../common/prepare_beds_and_contrast.R')
	#bedfiles are ready
	karyotype<-get.cytoband.ranges()
	#karyotype ready
	SNPchip.version=package.version('SNPchip')
	karyotype.with.methylation<-CountCoverageOfNoodles(karyotype,bedfilenames,DNAids)
	#it is folder with bed files
	save(file='karyotype.with.methylation.Rda',list=c('karyotype','karyotype.with.methylation','Clinical','clinFile','clinFileName','bedsinfolder','bed.used','tumors','normals','DNAids','SNPchip.version'))
}


karyotype.DM.loaded<-FALSE
# we can the whole thing to karyotype.with.methylation.Rda
if(file.exists('karyotype.DM.Rda'))
	if ('DM.cytobands' %in% load('karyotype.DM.Rda'))
			karyotype.DM.loaded<-TRUE
if(!karyotype.DM.loaded)
{
	tests.number<-dim(karyotype.with.methylation)[1]

	expected.w.statistic<-(sum(normals)*sum(tumors))/2

	wilcoxon.res<-apply(karyotype.with.methylation,1,function(row){
			meth.values<-jitter(row)
			w<-wilcox.test(meth.values[normals],meth.values[tumors])
			c(w$p.value,(w[['statistic']]<expected.w.statistic))
		})

	wilcoxon.p.values<-wilcoxon.res[1,]
	normals.are.less.methylated<-as.logical(wilcoxon.res[2,])

	wilcoxon.p.values.bonferroni<-p.adjust(wilcoxon.p.values,'bonferroni')
	wilcoxon.p.values.fdr<-p.adjust(wilcoxon.p.values,'fdr')

	DM.cytobands<-which(wilcoxon.p.values<=0.05)
	DM.cytobands.bonferroni<-which(wilcoxon.p.values.bonferroni<=0.05)
	DM.cytobands.fdr<-which(wilcoxon.p.values.fdr<=0.05)
	save(file='karyotype.DM.Rda',list=c('DM.cytobands','DM.cytobands.bonferroni','DM.cytobands.fdr'))
}
