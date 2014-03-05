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

data.loaded<-FALSE
# we can the whole thing to karyotype.with.methylation.Rda
if(file.exists('karyotype.with.methylation.Rda'))
	if ('karyotype.with.methylation' %in% load('karyotype.with.methylation.Rda'))
		if (class(karyotype.with.methylation)=='data.frame')
			data.loaded<-TRUE

if (!data.loaded)
{
	dataFolder<-'../../Data'

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

	DNAids=Clinical$'Tissue DNA HAND ID'
	#Clinical prepared.

	#it is folder with bed files
	peakbedsfolder<-paste0(dataFolder,'/PeakCalls/bedfiles/')

	beds<-list.files(peakbedsfolder)

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
	karyotype.with.methylation<-as(karyotype,"data.frame")

	bed_available<-logical(0)

	for (DNAid in DNAids)
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
		overlaps<-findOverlaps(karyotype,methylated.ranges)
		methylcoverage=integer(0)
		for(chr in names(karyotype))
		#cycle by chromosome
		{
			methylcoverage.this.chr<-sapply(1:length(karyotype[chr][[1]]),function(band){
				sum(width(methylated.ranges[chr][as.list(overlaps[[chr]])[[band]],]))
			})#list of methylated coverage per cytoband
			methylcoverage<-c(methylcoverage,methylcoverage.this.chr)
			#add to main list
			#we need dual cycle because of the findOverlap return structure
		}
		karyotype.with.methylation[[DNAid]]=methylcoverage
		bed_available<-c(bed_available,TRUE)
	}
	save(file='karyotype.with.methylation.Rda',list=c('karyotype.with.methylation','Clinical','bed_available','tumors','normals','DNAids'))
}

wilcoxon.p.values<-numeric(0)
normals.are.less.methylated<-logical(0)

expected.w.statistic<-(sum(normals[bed_available])*sum(tumors[bed_available]))/2

tests.number<-dim(karyotype.with.methylation)[1]

for (rown in 1:tests.number)
{
	meth.values<-as.numeric(karyotype.with.methylation[rown,][DNAids[bed_available]])
	w<-wilcox.test(meth.values[normals[bed_available]],meth.values[tumors[bed_available]])
	wilcoxon.p.values<-c(wilcoxon.p.values,w$p.value)
	normals.are.less.methylated<-c(normals.are.less.methylated,(w[['statistic']]<expected.w.statistic))
	#anova.result<-c(anova.result,anova(lm(meth.values~tumors))[1,'Pr(>F)'])
}

DM.cytobands<-which(wilcoxon.p.values<=0.05)
DM.cytobands.Bonferroni<-which(wilcoxon.p.values*tests.number<=0.05)

