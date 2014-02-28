library('rtracklayer')

dataFolder<-'../../Data'

#the file to be read
#clinFile <- paste0(dataFolder,'/Clinical/2013-11-25_RO1_Batch1_ClinicalData_DeID.xls.csv')
clinFile <- paste0(dataFolder,'/Clinical/2014-02-06_RO1_Batch1-4_ClinicalData_DeID.xls.csv')
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

length.sum<-integer(0)
score.sum<-numeric(0)
ls.prod.sum<-numeric(0)

#it is folder with bed files
peakbedsfolder<-paste0(dataFolder,'/PeakCalls/bedfiles/')

beds<-list.files(peakbedsfolder)

bed_available<-logical(0)

for (DNAid in Clinical$'Tissue DNA HAND ID')
{
		DNAidKey<-paste0(strsplit(DNAid,'_')[[1]],collapse='') #remove _ from key
		match<-grep(DNAidKey,beds)
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
		bed_available<-c(bed_available,TRUE)
}

Y.chrom.stats<-
	data.frame(
		code=Clinical$Code[bed_available],
		DNAid=Clinical$'Tissue DNA HAND ID'[bed_available],
		tumor=tumors[bed_available],
		male=Clinical$Gender[bed_available]=='1',
		'Peak length sum'=length.sum,
		'Peak score sum'=score.sum,
		'Product sum'=ls.prod.sum,
		stringsAsFactors = FALSE)

print(Y.chrom.stats)

print(wilcox.test(Y.chrom.stats$'Peak length sum'[Y.chrom.stats$male],Y.chrom.stats$'Peak length sum'[!Y.chrom.stats$male]))
print(wilcox.test(Y.chrom.stats$'Peak score sum'[Y.chrom.stats$male],Y.chrom.stats$'Peak score sum'[!Y.chrom.stats$male]))
print(wilcox.test(Y.chrom.stats$'Product sum'[Y.chrom.stats$male],Y.chrom.stats$'Product sum'[!Y.chrom.stats$male]))
