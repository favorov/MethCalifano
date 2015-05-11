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

if(!('all.for.huge.dasha.report.loaded' %in% ls()) || is.na(all.for.huge.dasha.report.loaded)) all.for.huge.dasha.report.loaded<-FALSE
#for quick-develop
if(!all.for.huge.dasha.report.loaded)
{
	message('loading or creating huge.report.frame')
	source('write.huge.report.prepare.R')
	message('loading methylation')
	load('noodles.C.Rda')
	message('loading normal reads')
	load('reads/noodles.C.7.spaghetti.normals.read.coverage.Rda')
	message('loading tumor reads')
	load('reads/noodles.C.7.spaghetti.tumors.read.coverage.Rda')
}

rows.no<-dim(noodles.C.methylation)[1]
report.interval<-1:rows.no


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
	report.framere<-huge.report.frame[fragment.range,]
	rownames(report.framere)=rownames(huge.report.frame)[fragment.range]
	meth.framere<-matrix(0,nrow=fragment.end-fragment.start+1,ncol=dim(noodles.C.methylation)[2])
	colnames(meth.framere)<-colnames(noodles.C.methylation)
	meth.framere[noodles.C.methylation[fragment.range,]>0]=1
	report.framere<-cbind(report.framere,meth.framere)
	readnormat<-as.matrix(noodles.C.7.spaghetti.normals.read.coverage[fragment.range,])
	colnames(readnormat)<-
		paste(colnames(noodles.C.7.spaghetti.normals.read.coverage),'.700.reads',sep='')
	report.framere<-cbind(report.framere,readnormat)
	readtummat<-as.matrix(noodles.C.7.spaghetti.tumors.read.coverage[fragment.range,])
	colnames(readtummat)<-
		paste(colnames(noodles.C.7.spaghetti.tumors.read.coverage),'.700.reads',sep='')
	report.framere<-cbind(report.framere,readtummat)
	colnames(report.framere)<-gsub(' ','',colnames(report.framere))
	write.table(report.framere,file=tsvfilename,sep='\t',quote=FALSE,row.names=TRUE,append=(fragment!=1),col.names=(fragment==1))
	rm(list=c('report.framere'))
	#header for first fragment
	#append for others
}
message('done')

