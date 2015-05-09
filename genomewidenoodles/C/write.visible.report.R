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

if (!suppressWarnings(require('xtable')))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite("xtable")
	library("xtable")
}

if (!suppressWarnings(require('Matrix')))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite("Matrix")
	library("Matrix")
}

if(!('all.the.all.loaded' %in% ls()) || is.na(all.the.all.loaded)) all.the.all.loaded<-FALSE
#for quick-develop
if(!all.the.all.loaded)
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
	all.the.all.loaded<-TRUE
}


generate.noodles.C.report<-function(report.set,#indices
												set.id, #variable part of the output file names
												no.html=FALSE) 
{
	report.noodles<-noodles.C[report.set,]
	report.fisher<-fisher.noodles.C.result[report.set,]
	rows.no<-length(report.set)

	tsvfilename=paste0("noodles.C.annotation.",set.id,".tsv")
	htmlfilename=paste0("noodles.C.annotation.",set.id,".html")


	#prepare dataframe
	message('init dataframe')
	report.frame<-data.frame('chr'=as.character(seqnames(report.noodles)),start=start(report.noodles),end=end(report.noodles),stringsAsFactors = FALSE)
		
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

	report.frame<-cbind(report.frame,elementMetadata(closest.genes)[,c('closest.TSS','pos','dir','dist')])

	message('done')

	message('Looking for overlapped genes')

	flanks<-7000

	ovelapped.genes<-genes.with.TSS.covered.by.interval(report.noodles,flanks=flanks)

	report.frame<-cbind(report.frame,elementMetadata(ovelapped.genes)[,c('overlapped.TSS','overlapped.pos','ovrl.dir')])

	message('done')

	
	message('Normal read stats')


	colnames(norm.read.stats.frame)<-c('norm.700.reads.min','norm.700.reads.25q','norm.700.reads.med','norm.700.reads.75q','norm.700.reads.max')

	report.frame<-cbind(report.frame,norm.read.stats.frame[report.set,])
	
	message('Tumor read stats')
	colnames(tumor.read.stats.frame)<-c('tumor.700.reads.min','tumor.700.reads.25q','tumor.700.reads.med','tumor.700.reads.75q','tumor.700.reads.max')
	
	report.frame<-cbind(report.frame,tumor.read.stats.frame[report.set,])

	message('done')
	#prepared
	
	save(file=paste0('noodles.C.annotation.',set.id,'.Rda'),list=c('report.frame'))
	
	write.table(report.frame,file=tsvfilename,sep='\t',row.names=FALSE,quote=FALSE)

	if(!no.html)
	{
		if(file.exists(htmlfilename)) {file.remove(htmlfilename)}

		print(xtable(report.frame,digits=c(0,0,0,0,8,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), display=c('d','s','d','d','g','f','f','f','f','f','s','s','s','s','s','d','s','d','s','s','s','d','d','d','d','d','d','d','d','d','d')), type="html", file=htmlfilename, include.rownames=FALSE)
#digits and display are to be +1 because of rows# that we do not print
	}
}

fish<-fisher.noodles.C.result$fisher.p.values

generate.noodles.C.report(which(p.adjust(fish,method='bonferroni')<=0.05),'bonf')
generate.noodles.C.report(which(p.adjust(fish,method='fdr')<=0.05),'fdr')

generate.noodles.C.report(which(fish<0.05),'uncorr',no.html=TRUE)

