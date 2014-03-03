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

DNAids=Clinical$'Tissue DNA HAND ID'
#Clinical prepared.
#reading karyotype

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

#print(karyotype)
#karyotype read

karyotype.with.counters<-karyotype

for (DNAid in DNAids)
{
	bedfilename<-paste0(dataFolder,'/PeakCalls/bedfiles/',DNAid,'_enrich.combined.removeDuplicates_vs_',DNAid,'_total.combined.removeDuplicates.50_peaks.bed')
	methylated.ranges<-as(import(bedfilename),"RangedData")
	overlaps<-findOverlaps(karyotype,methylated.ranges)
	methylcoverage=integer(0)
	for(chr in names(karyotype))
	{
		methylcoverage.this.chr<-sapply(1:length(karyotype[chr][[1]]),function(band){
			sum(width(methylated.ranges[chr][as.list(overlaps[[chr]])[[band]],]))
		})
		methylcoverage<-c(methylcoverage,methylcoverage.this.chr)
	}
	karyotype.with.counters[[as.character(DNAid)]]=methylcoverage
}


kdf<-as(karyotype.with.counters,'data.frame')

wilcoxon.result<-numeric(0)
anova.result<-numeric(0)

for (rown in 1:dim(kdf)[1])
{
	meth.values<-as.numeric(kdf[rown,paste0('X',DNAids)])
	wilcoxon.result<-c(wilcoxon.result,wilcox.test(meth.values[normals],meth.values[tumors])$p.value)
	anova.result<-c(anova.result,anova(lm(meth.values~tumors))[1,'Pr(>F)'])
}

karyotype$wilcoxon.result<-wilcoxon.result
karyotype$anova.result<-anova.result

DM.intervals<-which(karyotype$wilcoxon.result<=0.05 | karyotype$anova.result<=0.05)

print(karyotype[DM.intervals,])
print(karyotype.with.counters[DM.intervals,])

