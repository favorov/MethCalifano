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

if (!require('DASiR'))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite("DASiR")
}

vistaEnhancers.loaded<-FALSE
# we can the whole thing to karyotype.with.methylation.Rda

if(file.exists('vistaEnhancers.Rda'))
	if ('vistaEnhancers' %in% load('vistaEnhancers.Rda'))
		if (class(vistaEnhancers)=='GRanges')
			vistaEnhancers.loaded<-TRUE

if(!vistaEnhancers.loaded)
{
	setDasServer(server="http://genome.cse.ucsc.edu/cgi-bin/das/")
	source = "hg19"
	chrom.ranges<-getDasEntries(source,as.GRanges=TRUE)
	chrom.ranges<-chrom.ranges[nchar(as.character(seqnames(chrom.ranges))) < 3 & (as.character(seqnames(chrom.ranges))) != 'M' & (as.character(seqnames(chrom.ranges))) != 'Y']
	#vistaEnhancers does not work for Y
	#remove all pseudochromosomes: they have long names and remove MT
	#result is: chrom.ranges is all the chromosomes 1..22,X,Y with their length
	#print(chrom.ranges)
	vistaEnhancers_data<-getDasFeature(source,chrom.ranges,'vistaEnhancers')
	#vistaEnhancers is a DataFrame with all the Vista Enhacers bands enumerated
	chrs<-as.character(seqnames(chrom.ranges))[vistaEnhancers_data$segment.range]
	#segment.range in return is not the chrom name; it is index of the chr's record in chrom.range
	vistaEnhancers<-GRanges(
		seqnames=paste0('chr',chrs), 
		ranges=IRanges
		(
			start=as.numeric(as.character(vistaEnhancers_data$start)),
			end=as.numeric(as.character(vistaEnhancers_data$end))
		),
		seqinfo=nucl.chromosomes.hg19(),
		id=vistaEnhancers_data$label,
		score=as.numeric(as.character(vistaEnhancers_data$score)),
		type=as.character(vistaEnhancers_data$type)
	)
	vista.enhancers.timestamp<-Sys.time()
  #log of conversion. accomodate and test before use
  #ve=as(vistaEnhancers,'GRanges')
  #elementMetadata(ve)$id=as.character(elementMetadata(ve)$id)
  #elementMetadata(ve)$score=as.integer(as.character(elementMetadata(ve)$score))
  #vistaEnhancers<-ve
  #up to here
	save(file='vistaEnhancers.Rda',list=c('vistaEnhancers','vista.enhancers.timestamp'))
}
