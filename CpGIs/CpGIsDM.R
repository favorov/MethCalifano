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

CpGIs.with.methylation.loaded<-FALSE
# we can the whole thing to CpGIs.with.methylation.Rda
if(file.exists('CpGIs.with.methylation.Rda'))
	if ('CpGIs.with.methylation' %in% load('CpGIs.with.methylation.Rda'))
		if (class(CpGIs.with.methylation)=='data.frame')
			CpGIs.with.methylation.loaded-TRUE

if (!CpGIs.with.methylation.loaded)
{
	dataFolder<-'../../Data'

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

	DNAids=Clinical$'Tissue DNA HAND ID'
	#Clinical prepared.

	#it is folder with bed files
	peakbedsfolder<-paste0(dataFolder,'/PeakCalls/bedfiles/')

	beds<-list.files(peakbedsfolder)

	#reading islands 
	# we can the whole thing to CpGIs.with.methylation.Rda
	source('load_or_read_and_save_CpGIs.R')
	#islands are read fro DAS or loaded
	CpGIs.with.methylation<-as(CpGIs,"data.frame")

	bed_available<-logical(0)
	bed_used<-rep(FALSE,length(beds))

	for (DNAid in DNAids)
	{
		DNAidKey<-strsplit(DNAid,',')[[1]][1]	#remove all after ,	
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
		overlaps<-findOverlaps(CpGIs,methylated.ranges)
		methylcoverage=integer(0)
		for(chr in names(CpGIs))
		#cycle by chromosome
		{
			methylcoverage.this.chr<-sapply(1:length(CpGIs[chr][[1]]),function(band){
				sum(width(methylated.ranges[chr][as.list(overlaps[[chr]])[[band]],]))
			})#list of methylated coverage per isaland 
			methylcoverage<-c(methylcoverage,methylcoverage.this.chr)
			#add to main list
			#we need dual cycle because of the findOverlap return structure
		}
		CpGIs.with.methylation[[DNAid]]=methylcoverage
		bed_used[match[1]]<-TRUE
		bed_available<-c(bed_available,TRUE)
	}
	save(file='CpGIs.with.methylation.Rda',list=c('CpGIs.with.methylation','Clinical','clinFile','beds','bed_available','bed_used','tumors','normals','DNAids'))
}

wilcoxon.p.values<-numeric(0)
normals.are.less.methylated<-logical(0)

expected.w.statistic<-(sum(normals[bed_available])*sum(tumors[bed_available]))/2

tests.number<-dim(CpGIs.with.methylation)[1]

for (rown in 1:tests.number)
{
	meth.values<-as.numeric(CpGIs.with.methylation[rown,][DNAids[bed_available]])
	w<-wilcox.test(meth.values[normals[bed_available]],meth.values[tumors[bed_available]])
	wilcoxon.p.values<-c(wilcoxon.p.values,w$p.value)
	normals.are.less.methylated<-c(normals.are.less.methylated,(w[['statistic']]<expected.w.statistic))
	#anova.result<-c(anova.result,anova(lm(meth.values~tumors))[1,'Pr(>F)'])
}

DM.islands<-which(wilcoxon.p.values<=0.05)
DM.islands.Bonferroni<-which(wilcoxon.p.values*tests.number<=0.05)

