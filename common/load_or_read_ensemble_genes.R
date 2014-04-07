#fetches ensemble genes from UCSC DAS server
#saves the information to ../common/karyotype.Rda
#if the file exists, read it after some checks
#From UCSC FAQ:
##Linking gene name with accession number	
#Question: 
#"I have the accession number for a gene and would like to link it to the gene name. Is there a table that shows both pieces of information?"
#Response:
#If you are looking at the RefSeq Genes, the refFlat table contains both the gene name (usually a HUGO Gene Nomenclature Committee ID) and its accession number. For the Known Genes, use the kgAlias table.

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

ensGenes.loaded<-FALSE
# we can the whole thing to ensGenes.with.methylation.Rda

if(file.exists('../common/ensGenes.Rda'))
	if ('ensGenes' %in% load('../common/ensGenes.Rda'))
		if (class(ensGenes)=='RangedData')
			ensGenes.loaded<-TRUE

if(!ensGenes.loaded)
{
	setDasServer(server="http://genome.cse.ucsc.edu/cgi-bin/das/")
	source = "hg19"
	chrom.ranges<-getDasEntries(source,as.GRanges=TRUE)
	chrom.ranges<-chrom.ranges[nchar(as.character(seqnames(chrom.ranges))) < 3 & (as.character(seqnames(chrom.ranges))) != 'M']
	#remove all pseudochromosomes: they have long names and remove MT
	#result is: chrom.ranges is all the chromosomes 1..22,X,Y with their length
	#print(chrom.ranges)
	ensGenes_data<-getDasFeature('hg19',chrom.ranges,'ensGene')
	chrs<-as.character(seqnames(chrom.ranges))[ensGenes_data$segment.range]
	#segment.range in return is not the chrom name; it is index of the chr's record in chrom.range
	ensGenes<-RangedData(
		space=paste0('chr',chrs), 
		ranges=IRanges
		(
			start=as.numeric(as.character(ensGenes_data$start)),
			end=as.numeric(as.character(ensGenes_data$end))
		),
		id=ensGenes_data$id,
		label=ensGenes_data$label,
		orientation=ensGenes_data$orientation
	)
	save(file='../common/ensGenes.Rda',list=c('ensGenes'))
}

