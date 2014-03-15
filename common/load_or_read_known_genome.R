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

setDasServer(server="http://genome.cse.ucsc.edu/cgi-bin/das/")
source = "hg19"
gr<-GRanges(seqnames =
	Rle(c("chr1"), c(1)),
	ranges = IRanges(204131200,204135400))
regs<-getDasFeature('hg19',gr,'knownGene')
print (regs)

