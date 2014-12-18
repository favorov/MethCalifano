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

#if (!suppressWarnings(require('caTools')))
#{
#	source("http://bioconductor.org/biocLite.R")
#	biocLite("caTools")
#	library("caTools")
#}

if(!('all.the.all.loaded' %in% ls()) || is.na(all.the.all.loaded)) all.the.all.loaded<-FALSE
#for quick-develop
if(!all.the.all.loaded)
{
	message('loading..')
	load('noodles.C.Rda')
	load('noodles.C.fisher.results.Rda')
	load('reads/noodles.C.7.spaghetti.normals.read.coverage.Rda')
	load('../../CytoBands/cytobands.DM.Rda')
	load('../../CpGIs/CpGIs.Rda')
	load('../../CpGIs/CpGIs.DM.indices.Rda')
	all.the.all.loaded<-TRUE
}

rows.no<-dim(fisher.noodles.C.result)[1]
report.interval<-1:rows.no


huge.loaded<-FALSE
if(file.exists('huge.Rda'))
	if ('report.frame' %in% load('huge.Rda'))
		if (class(report.frame)=='data.frame')
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

	message('Normal read stats')
	
	#spaghetti.size.in.noodles<-7
	
	#spaghetti.C.normals.read.coverage<-
	#	spaghetti.size.in.noodles*
	#	caTools::runmean(noodles.C.normals.read.coverage,spaghetti.size.in.noodles,alg='fast')
	#running mean*window.size is running sum
	#S4Vectors::runmean tries to shade the caTools::runmean

	norm.read.stats.frame<-t(apply(noodles.C.7.spaghetti.normals.read.coverage[report.interval,],1,quantile))

	colnames(norm.read.stats.frame)<-c('norm.700.reads.min','norm.700.reads.25q','norm.700.reads.med','norm.700.reads.75q','norm.700.reads.max')

	report.frame<-cbind(report.frame,norm.read.stats.frame)

	message('done')
	#prepared

	save(file='huge.Rda',list=c('report.frame'))
}

message('writing')

tsvfilename="noodles.C.complete.annotaion.tsv"
fragments.to.out<-50
step<-rows.no %/% fragments.to.out + ifelse(rows.no %% fragments.to.out > 0,1,0) #if remainder is zero, / is ok 

for(fragment in 1:fragments.to.out)
{	
	message(paste0('Fragment ',fragment, ' of ',fragments.to.out))
	fragment.start <- 1 + step*(fragment-1)
	fragment.end <- min(fragment.start+step-1,rows.no)
	fragment.range<-fragment.start:fragment.end
	report.framere<-report.frame[fragment.range,]
	rownames(report.framere)=rownames(report.frame)[fragment.range]
	meth.framere<-matrix(0,nrow=fragment.end-fragment.start+1,ncol=dim(noodles.C.methylation)[2])
	colnames(meth.framere)<-colnames(noodles.C.methylation)
	meth.framere[noodles.C.methylation[fragment.range,]>0]=1
	report.framere<-cbind(report.framere,meth.framere)
	readmat<-noodles.C.7.spaghetti.normals.read.coverage[fragment.range]
	colnames(readmat)<-
		paste(colnames(noodles.C.7.spaghetti.normals.read.coverage),'.700.reads',sep='')
	report.framere<-cbind(report.framere,readmat)
	colnames(report.framere)<-gsub(' ','',colnames(report.framere))
	write.table(report.framere,file=tsvfilename,sep='\t',quote=FALSE,row.names=TRUE,append=(fragment!=1),col.names=(fragment==1))
	#header for first fragment
	#append for others
}

