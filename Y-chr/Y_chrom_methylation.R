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
	length.sum<-integer(0)
	score.sum<-numeric(0)
	ls.prod.sum<-numeric(0)

	#it is folder with bed files
	peakbedsfolder<-paste0(dataFolder,'/PeakCalls/bedfiles/')

	beds<-list.files(peakbedsfolder)

	bed_available<-logical(0)
	bed_used<-rep(FALSE,length(beds))

	for (DNAid in DNAids)
	{
			DNAidKey<-strsplit(DNAid,',')[[1]][1]		#remove all after ,	
			match<-grep(DNAidKey,beds)
			if (!length(match)) 
			{
				DNAidKey<-paste0(strsplit(DNAid,'_')[[1]],collapse='') 
				#remove _ from key; sometimes, it help
				match<-grep(DNAidKey,beds)
			}
			if (!length(match)) 
			{
				bed_available<-c(bed_available,FALSE)
				next
			}
			if (length(match)>1) stop(paste0("More than one match of DNAid ",DNAid," amonng the bed file names.\n"));
			bedfilename<-paste0(peakbedsfolder,beds[match[1]]);
			methylated.ranges<-as(import(bedfilename),"RangedData")
			length.sum<-c(length.sum,sum(width(methylated.ranges['chrY'])))
			score.sum<-c(score.sum,sum(methylated.ranges['chrY']$score))
			ls.prod.sum<-c(ls.prod.sum,sum(width(methylated.ranges['chrY'])*methylated.ranges['chrY']$score))
			bed_used[match[1]]<-TRUE
			bed_available<-c(bed_available,TRUE)
	}

	Y.chrom.stats<-
		data.frame(
			code=Clinical$Code[bed_available],
			DNAid=DNAids[bed_available],
			tumor=tumors[bed_available],
			male=Clinical$Gender[bed_available]=='1',
			'Peak length sum'=length.sum,
			'Peak score sum'=score.sum,
			'Product sum'=ls.prod.sum,
			stringsAsFactors = FALSE)
	save(file='Y.chrom.stats.Rda',list=c('Y.chrom.stats','clinFile','Clinical','beds','bed_available','bed_used','DNAids'))

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

