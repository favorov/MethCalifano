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

if (!suppressWarnings(require('Matrix')))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite("Matrix")
	library("Matrix")
}

if (!suppressWarnings(require('data.table')))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite("data.table")
	library("data.table")
}

if(!('all.for.visible.dasha.report.loaded' %in% ls()) || is.na(all.for.visible.dasha.report.loaded)) all.for.visible.dasha.report.loaded<-FALSE
#for quick-develop
if(!all.for.visible.dasha.report.loaded)
{
	message('loading..')
	load('noodles.C.Rda')
	load('noodles.C.fisher.results.Rda')
	load('reads/noodles.C.7.spaghetti.normals.read.quantiles.Rda')
	load('reads/noodles.C.7.spaghetti.tumors.read.quantiles.Rda')
	load('xeno.C.methylation.Rda')
	load('../../CytoBands/cytobands.DM.Rda')
	load('../../CpGIs/CpGIs.Rda')
	load('../../CpGIs/CpGIs.DM.indices.Rda')
	colnames(norm.read.stats.frame)<-c('norm.700.reads.min','norm.700.reads.25q','norm.700.reads.med','norm.700.reads.75q','norm.700.reads.max')
	colnames(tumor.read.stats.frame)<-c('tumor.700.reads.min','tumor.700.reads.25q','tumor.700.reads.med','tumor.700.reads.75q','tumor.700.reads.max')
	all.for.visible.dasha.report.loaded<-TRUE
}

rows.no<-dim(fisher.noodles.C.result)[1]
report.interval<-1:rows.no


huge.loaded<-FALSE
if(file.exists('huge.report.frame.Rda'))
	if ('huge.report.frame' %in% load('huge.report.frame.Rda'))
		if (class(huge.report.frame)=='data.frame')
			huge.loaded<-TRUE

if(!huge.loaded)
{
	report.noodles<-noodles.C[report.interval,]
	report.fisher<-fisher.noodles.C.result[report.interval,]
	#actually, it is to develop for little tests.no

	#prepare dataframe
	message('init dataframe')
	report.frame<-data.frame('chr'=as.character(seqnames(report.noodles)),start=start(report.noodles),end=end(report.noodles),stringsAsFactors = FALSE)

	rownames(report.frame)<-paste(report.frame$chr,":",report.frame$start,'-',report.frame$end,sep='')

	message('adding Fisher')
	report.frame<-cbind(report.frame,
			'fisher.p.value'=report.fisher$fisher.p.values,
			'tmr.ratio'=report.fisher$meth.in.tumors.ratio,
			'nor.ratio'=report.fisher$meth.in.normals.ratio,
			'OR'=report.fisher$OR,
			'CI_95_L'=report.fisher$CI_95_L,
			'CI_95_H'=report.fisher$CI_95_H
		)

	message('Mapping to karyotype...')

	cb<-integer(rows.no)

	noodles.to.karyotype<-findOverlaps(report.noodles,cytobands,type="within")

	cb[queryHits(noodles.to.karyotype)]=subjectHits(noodles.to.karyotype)

	cb[cb==0]=NA

	message('done')

	report.frame<-cbind(report.frame,'cytoband'=cytobands$'name'[cb],'DM.band?'=cytobands.DM.statistics$'wilcoxon.p.values'[cb]<0.05,stringsAsFactors = FALSE)
	#prepared

	message('Mapping to cpg islands...')

	noodles.to.cpgi<-findOverlaps(report.noodles,CpGIs,type="within")

	ci<-integer(rows.no)

	ci[queryHits(noodles.to.cpgi)]=subjectHits(noodles.to.cpgi)

	ci[ci==0]=NA

	message('done')

	report.frame<-cbind(report.frame,'CpGi'=CpGIs$'id'[ci],'DM.island?'=ifelse(is.na(ci),NA,as.logical(ci %in% DM.CpGIslands)),stringsAsFactors = FALSE)

	#report.frame$'CpGi'<-substr(report.frame$'CpGi',6,1000) # 1000 'any'; we strip first 'CpGi: ' from the id

	message('Looking for closest genes')
	closest.genes<-closest.gene.start.by.interval(report.noodles)

	message('done')


	message('combining')

	report.frame<-cbind(report.frame,elementMetadata(closest.genes)[,c('closest.TSS','pos','dir','dist')])
	message('done\n')

	message('done')

	message('Xeno and cell lines')

	xeno<-data.frame(matrix(ncol=4,ifelse(xeno.C.methylation[report.interval,]>0,1,0)))

	colnames(xeno)<-colnames(xeno.C.methylation)

	report.frame<-cbind(report.frame,xeno)

	message('done')

	#prepared

	huge.report.frame<-data.frame(report.frame)

	message('writing')
	
	save(file='huge.report.frame.Rda',list=c('huge.report.frame'))
	
	message('done')
}


