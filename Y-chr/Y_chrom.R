library('rtracklayer')

dataFolder<-'../../Data'

clinFile <- paste0(dataFolder,'/Clinical/2013-11-25_RO1_Batch1_ClinicalData_DeID.xls.csv')
Clinical <- read.csv(clinFile,stringsAsFactors=F)

# fix up the column names to account for the fact that there are two header rows
colnames(Clinical)[1:3] <- Clinical[1,1:3]
Clinical <- Clinical[2:nrow(Clinical),]

#remove all service rows

Clinical<-Clinical[!is.na(as.integer(Clinical[,3])),]

# determine which are tumor and which are normal samples
tumors <- as.integer(Clinical$Code) < 99
normals <- as.integer(Clinical$Code) > 99

Y.chrom.stats<-data.frame(code=Clinical$Code,DNAid=Clinical$'Tissue DNA HAND ID',tumor=tumors,male=Clinical$Gender=='1',stringsAsFactors = FALSE)

length.sum<-integer(0)
score.sum<-numeric(0)
ls.prod.sum<-numeric(0)

for (DNAid in Y.chrom.stats$DNAid)
{
	bedfilename<-paste0(dataFolder,'/PeakCalls/bedfiles/',DNAid,'_enrich.combined.removeDuplicates_vs_',DNAid,'_total.combined.removeDuplicates.50_peaks.bed')
	methylated.ranges<-as(import(bedfilename),"RangedData")
	length.sum<-c(length.sum,sum(width(methylated.ranges['chrY'])))
	score.sum<-c(score.sum,sum(methylated.ranges['chrY']$score))
	ls.prod.sum<-c(ls.prod.sum,sum(width(methylated.ranges['chrY'])*methylated.ranges['chrY']$score))
}

Y.chrom.stats<-cbind(Y.chrom.stats,'Peak length sum'=length.sum,'Peak score sum'=score.sum,'Product sum'=ls.prod.sum)

print(Y.chrom.stats)

print(wilcox.test(Y.chrom.stats$'Peak length sum'[Y.chrom.stats$male],Y.chrom.stats$'Peak length sum'[!Y.chrom.stats$male]))
print(wilcox.test(Y.chrom.stats$'Peak score sum'[Y.chrom.stats$male],Y.chrom.stats$'Peak score sum'[!Y.chrom.stats$male]))
print(wilcox.test(Y.chrom.stats$'Product sum'[Y.chrom.stats$male],Y.chrom.stats$'Product sum'[!Y.chrom.stats$male]))
