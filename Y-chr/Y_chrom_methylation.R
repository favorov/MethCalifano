if (!require('rtracklayer'))
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
	dataFolder<-'../../Data'

	#the file to be read
	#clinFile <- paste0(dataFolder,'/Clinical/2013-11-25_RO1_Batch1_ClinicalData_DeID.xls.csv')
	clinFile <- paste0(dataFolder,'/Clinical/2014-02-06_RO1_Batch1-4_ClinicalData_DeID.csv')
	Clinical <- read.csv(clinFile,stringsAsFactors=F)

	#fix up the column names to account for the fact that there are two header rows
	#the names are absent for first three rows, they were in a row below 
	colnames(Clinical)[1:3] <- Clinical[1,1:3]
	Clinical <- Clinical[2:nrow(Clinical),]

	#remove all service rows (where the TissueRNAID is not integer id)
	Clinical<-Clinical[!is.na(as.integer(Clinical[,3])),]

	#rownames are Codes
	rownames(Clinical)<-Clinical$Codes

	# determine which are tumor and which are normal samples
	tumors <- as.integer(Clinical$Code) < 99
	normals <- as.integer(Clinical$Code) > 99

	DNAids<-Clinical$'Tissue DNA HAND ID'

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

