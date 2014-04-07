#fetches known CpGIs from UCSC DAS server
#saves the information to ../common/CpGIs.Rda
#if the file exists, read it after some checks

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

CpGIs.loaded<-FALSE
# we can the whole thing to karyotype.with.methylation.Rda

if(file.exists('../common/CpGIs.Rda'))
	if ('CpGIs' %in% load('../common/CpGIs.Rda'))
		if (class(CpGIs)=='RangedData')
			CpGIs.loaded<-TRUE

if(!CpGIs.loaded)
{

	setDasServer(server="http://genome.cse.ucsc.edu/cgi-bin/das/")
	source = "hg19"
	chrom.ranges<-getDasEntries(source,as.GRanges=TRUE)
	chrom.ranges<-chrom.ranges[nchar(as.character(seqnames(chrom.ranges))) < 3 & (as.character(seqnames(chrom.ranges))) != 'M']
	#remove all pseudochromosomes: they have long names and remove MT
	#result is: chrom.ranges is all the chromosomes 1..22,X,Y with their length
	#print(chrom.ranges)
	CpGI_data<-getDasFeature(source,chrom.ranges,'cpgIslandExt')
	#CpGI_data is a DataFrame with all the CpGIs bands enumerated
	chrs<-as.character(seqnames(chrom.ranges))[CpGI_data$segment.range]
	#the same trick
	CpGIs<-RangedData(
		space=paste0('chr',chrs), 
		ranges=IRanges
		(
			start=as.numeric(as.character(CpGI_data$start)),
			end=as.numeric(as.character(CpGI_data$end))
		),
		id=CpGI_data$id
	)
	save(file='../common/CpGIs.Rda',list=c('CpGIs'))
}
#print(CpGIs)

