if (!suppressWarnings(require('rtracklayer')))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite("rtracklayer")
}

Y.chrom.stats.loaded<-FALSE
# we can the whole thing to karyotype.with.methylation.Rda
if(file.exists('Y.chrom.stats.Rda'))
	if ('Y.chrom.stats' %in% load('Y.chrom.stats.Rda'))
		if (class(Y.chrom.stats)=='data.frame')
			Y.chrom.stats.loaded<-TRUE
if (!Y.chrom.stats.loaded)
{
	source('../common/read_clinical.R')
	source('../common/prepare_beds_and_conrast.R')

	length.sum<-integer(0)
	score.sum<-numeric(0)
	ls.prod.sum<-numeric(0)

	for (bedfilename in bedfilenames)
	{
			methylated.ranges<-as(import(bedfilename),"RangedData")
			length.sum<-c(length.sum,sum(width(methylated.ranges['chrY'])))
			score.sum<-c(score.sum,sum(methylated.ranges['chrY']$score))
			ls.prod.sum<-c(ls.prod.sum,sum(width(methylated.ranges['chrY'])*methylated.ranges['chrY']$score))
	}

	Y.chrom.stats<-
		data.frame(
			code=codes,
			DNAid=DNAids,
			tumor=contrast,
			male=Clinical$Gender[clinical.row.used]=='1',
			'Peak length sum'=length.sum,
			'Peak score sum'=score.sum,
			'Product sum'=ls.prod.sum,
			stringsAsFactors = FALSE)
	save(file='Y.chrom.stats.Rda',list=c('Y.chrom.stats','clinFile','Clinical','bedsinfolder','bedfilenames','clinical.row.used','bed.used','DNAids'))
}

print(Y.chrom.stats)

#to compare with pure 1's and 0's: 
one_m<-rep(1,sum(Y.chrom.stats$male))
one_m_j<-jitter(one_m)
zero_f<-rep(0,sum(!Y.chrom.stats$male))
zero_f_j<-jitter(zero_f)

Y.chrom.stats$Peak.length.sum<-jitter(Y.chrom.stats$Peak.length.sum)
Y.chrom.stats$Peak.score.sum<-jitter(Y.chrom.stats$Peak.score.sum)
Y.chrom.stats$Product.sum<-jitter(Y.chrom.stats$Product.sum)

