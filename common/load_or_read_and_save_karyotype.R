#fetches cytobands from UCSC DAS server
#saves the information to ../common/karyotype.Rda
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

karyotype.loaded<-FALSE
# we can the whole thing to karyotype.with.methylation.Rda

if(file.exists('../common/karyotype.Rda'))
	if ('karyotype' %in% load('../common/karyotype.Rda'))
		if (class(karyotype)=='RangedData')
			karyotype.loaded<-TRUE

if(!karyotype.loaded)
{
	setDasServer(server="http://genome.cse.ucsc.edu/cgi-bin/das/")
	source = "hg19"
	chrom.ranges<-getDasEntries(source,as.GRanges=TRUE)
	chrom.ranges<-chrom.ranges[nchar(as.character(seqnames(chrom.ranges))) < 3 & (as.character(seqnames(chrom.ranges))) != 'M']
	#remove all pseudochromosomes: they have long names and remove MT
	#result is: chrom.ranges is all the chromosomes 1..22,X,Y with their length
	#print(chrom.ranges)
	karyotype_data<-getDasFeature(source,chrom.ranges,'cytoBand')
	#karyotype is a DataFrame with all the chromosome karyotype bands enumerated
	chrs<-as.character(seqnames(chrom.ranges))[karyotype_data$segment.range]
	#segment.range in return is not the chrom name; it is index of the chr's record in chrom.range
	karyotype<-RangedData(
		space=paste0('chr',chrs), 
		ranges=IRanges
		(
			start=as.numeric(as.character(karyotype_data$start)),
			end=as.numeric(as.character(karyotype_data$end))
		),
		id=karyotype_data$label
	)
	save(file='../common/karyotype.Rda',list=c('karyotype'))
}
#print(karyotype)

