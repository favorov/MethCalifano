if (!require('DASiR'))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite("DASiR")
}

setDasServer(server="http://genome.cse.ucsc.edu/cgi-bin/das/")
source = "hg19"
chrom.ranges<-getDasEntries(source,as.GRanges=TRUE)
chrom.ranges<-chrom.ranges[nchar(as.character(seqnames(chrom.ranges))) < 3 & (as.character(seqnames(chrom.ranges))) != 'M']
#remove all pseudochromosomes: they have long names and remove MT
#result is: chrom.ranges is all the chromosomes 1..22,X,Y with their length
#print(chrom.ranges)
cpgi_data<-getDasFeature(source,chrom.ranges,'cpgIslandExt')
#cpgi_data is a DataFrame with all the cpgis bands enumerated
chrs<-as.character(seqnames(chrom.ranges))[cpgi_data$segment.range]
#the same trick
cpgis<-RangedData(
	space=paste0('chr',chrs), 
	ranges=IRanges
	(
		start=as.numeric(as.character(cpgi_data$start)),
		end=as.numeric(as.character(cpgi_data$end))
	),
	id=cpgi_data$id
)
#now, we want to know how much cpgis each cytoband covers
print(cpgis)
