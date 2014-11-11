#fetches chrom ranges from UCSC DAS server
#saves the information to ../common/chromRanges.Rda
#if the file exists, read it after some checks
#From UCSC FAQ:
##Linking gene name with accession number	
#Question: 
#"I have the accession number for a gene and would like to link it to the gene name. Is there a table that shows both pieces of information?"
#Response:
#If you are looking at the RefSeq Genes, the refFlat table contains both the gene name (usually a HUGO Gene Nomenclature Committee ID) and its accession number. For the Known Genes, use the kgAlias table.

#if (!require('rtracklayer'))
#{
#	source("http://bioconductor.org/biocLite.R")
#	biocLite("rtracklayer")
#}

.Deprecated(new='get.cytodand.ranges',old='load_or_read_chrom_ranges.R')

if (!require('DASiR'))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite("DASiR")
}

chrom.ranges.loaded<-FALSE
# we can the whole thing to refseqGenes.with.methylation.Rda

if(file.exists('../common/chromRanges.Rda'))
	if ('chrom.ranges' %in% load('../common/chromRanges.Rda'))
			chrom.ranges.loaded<-TRUE

if(!chrom.ranges.loaded)
{
	setDasServer(server="http://genome.cse.ucsc.edu/cgi-bin/das/")
	source <- "hg19"
	chrom.ranges<-getDasEntries(source,as.GRanges=TRUE)
	chrom.ranges<-chrom.ranges[nchar(as.character(seqnames(chrom.ranges))) < 3 & (as.character(seqnames(chrom.ranges))) != 'M']
	#we removed the strange and test chromosomes and MT
	chrom.length <- end(ranges(chrom.ranges))
	names(chrom.length) <- paste0('chr', seqnames(chrom.ranges)) 
	save(file='../common/chromRanges.Rda',list=c('chrom.ranges','chrom.length'))
}
