#fetches cytobands from UCSC DAS server
#saves the information to ../common/karyotype.Rda
#if the file exists, read it after some checks

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

if (!require('rtracklayer'))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite("rtracklayer")
}



superEnhancers.loaded<-FALSE
# we can the whole thing to karyotype.with.methylation.Rda

if(file.exists('superEnhancers.Rda'))
	if ('superEnhancers' %in% load('superEnhancers.Rda'))
		if (class(superEnhancers)=='GRanges')
			superEnhancers.loaded<-TRUE

if(!superEnhancers.loaded)
{
	superDatafolder<-('../../Data/superEnhancers/')
	super.enhancers.timestamp<-Sys.time()
	dirbed<-dir(superDatafolder)
	dirbed<-dirbed[grep('bed',dirbed)]
	dirbed<-paste(superDatafolder,dirbed,sep='')
	superEnhancers<-GRanges(seqinfo=nucl.chromosomes.hg19(),type=character(0),is.super=logical(0))
	for (supenhbed in dirbed)
	{
		twobed<-import(supenhbed,seqinfo=nucl.chromosomes.hg19())
		bed<-twobed[[1]]
		bed$name<-NULL
		bed$type<-strsplit(names(twobed)[1],split=' ')[[1]][3]
		shortlist<-twobed[[2]]
		bed$short.list<-overlapsAny(bed,shortlist,type='within')
  	superEnhancers=c(superEnhancers,bed)
	}
	superEnhancers=superEnhancers[order(superEnhancers)]
	save(file='superEnhancers.Rda',list=c('superEnhancers','super.enhancers.timestamp'))
}
