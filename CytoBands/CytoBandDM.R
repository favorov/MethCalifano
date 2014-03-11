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

karyotype.with.methylation.loaded<-FALSE
# we can the whole thing to karyotype.with.methylation.Rda
if(file.exists('karyotype.with.methylation.Rda'))
	if ('karyotype.with.methylation' %in% load('karyotype.with.methylation.Rda'))
		if (class(karyotype.with.methylation)=='data.frame')
			karyotype.with.methylation.loaded<-TRUE

if (!karyotype.with.methylation.loaded)
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

	#reading karyotype
	# we can the whole thing to karyotype.with.methylation.Rda
	source('load_or_read_and_save_karyotype.R')
	#print(karyotype)
	#karyotype read fro DAS or loaded
	karyotype.with.methylation<-as(karyotype,"data.frame")

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
		bed_used[match[1]]<-TRUE
		bed_available<-c(bed_available,TRUE)
	}
	save(file='karyotype.with.methylation.Rda',list=c('karyotype.with.methylation','Clinical','clinFile','beds','bed_available','bed_used','tumors','normals','DNAids'))
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

wilcoxon.p.values.bonferroni<-p.adjust(wilcoxon.p.values,'bonferroni')
wilcoxon.p.values.fdr<-p.adjust(wilcoxon.p.values,'fdr')

DM.cytobands<-which(wilcoxon.p.values<=0.05)
DM.cytobands.bonferroni<-which(wilcoxon.p.values.bonferroni<=0.05)
DM.cytobands.fdr<-which(wilcoxon.p.values.fdr<=0.05)

