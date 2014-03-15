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
			CpGIs.with.methylation.loaded<-TRUE

if (!CpGIs.with.methylation.loaded)
{
	source('../common/read_clinical.R')
	#Clinical prepared.

	#it is folder with bed files
	peakbedsfolder<-paste0(dataFolder,'/PeakCalls/bedfiles/')

	beds<-list.files(peakbedsfolder)

	#reading islands 
	# we can the whole thing to CpGIs.with.methylation.Rda
	source('../common/load_or_read_and_save_CpGIs.R')
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
			})
			#list of methylated coverage per island
			#the main idead is that as.list(hitsobject) convert is to list of vectors, 
			#list is addressd by query oblects and give ref object lists (possibly, empty)
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

CpGIs.wilcoxon.data.loaded<-FALSE
# we can the whole thing to CpGIs.with.methylation.Rda
if(file.exists('CpGIs.wilcoxon.data.Rda'))
	if ('DM.CpGIslands.Bonferroni' %in% load('CpGIs.wilcoxon.data.Rda'))
			CpGIs.wilcoxon.data.loaded<-TRUE
if(!CpGIs.wilcoxon.data.loaded)
{
	wilcoxon.p.values<-numeric(0)

	normals.are.less.methylated<-logical(0)

	expected.w.statistic<-(sum(normals[bed_available])*sum(tumors[bed_available]))/2

	tests.number<-dim(CpGIs.with.methylation)[1]

	for (rown in 1:tests.number)
	{
		meth.values<-jitter(as.numeric(CpGIs.with.methylation[rown,][DNAids[bed_available]]))
		w<-wilcox.test(meth.values[normals[bed_available]],meth.values[tumors[bed_available]])
		wilcoxon.p.values<-c(wilcoxon.p.values,w$p.value)
		normals.are.less.methylated<-c(normals.are.less.methylated,(w[['statistic']]<expected.w.statistic))
		#anova.result<-c(anova.result,anova(lm(meth.values~tumors))[1,'Pr(>F)'])
	}

	DM.CpGIslands<-which(wilcoxon.p.values<=0.05)
	DM.CpGIslands.Bonferroni<-which(wilcoxon.p.values*tests.number<=0.05)
	save(file='CpGIs.wilcoxon.data.Rda',list=c('wilcoxon.p.values','normals.are.less.methylated','tests.number','DM.CpGIslands','DM.CpGIslands.Bonferroni'))
}

#here, we form output statictics
columns<-c('id','space','start','end')
CpGIs.stat<-cbind(CpGIs.with.methylation[,columns],'wilcoxon.p.value'=wilcoxon.p.values,'is.hyper'=normals.are.less.methylated)
DM.CpGIs.stat<-CpGIs.stat[DM.CpGIslands.Bonferroni,]

#we want to put each diffmet CpGi to a cytoband
source('load_or_read_and_save_karyotype.R')
DM.CpGIs.Ranges<-as(DM.CpGIs.stat[,columns],'RangedData')
CpGIs.to.karyotype<-findOverlaps(DM.CpGIs.Ranges,karyotype,type="within")
cytobands.of.DM.cpgis=character(0)
for (chr in names(CpGIs.to.karyotype))
{
  #chromosome cycle
	len<-length(DM.CpGIs.Ranges[chr][[1]])
	if (len==0) next
	#we did not do it in 'big' cycles over
	#overlaps, because their queries were like cytobands, etc - no empty chromosomes
	cytoband.numbers.this.chr<-sapply(1:len,function(island_no){
		my_cytoz<-as.list(CpGIs.to.karyotype[[chr]])[[island_no]]
		if (!length(my_cytoz)) {NA} else {my_cytoz[1]}
	})
	cytobands.of.DM.cpgis<-c(cytobands.of.DM.cpgis,as.character(karyotype[chr][[1]][cytoband.numbers.this.chr]))
}

DM.CpGIs.stat<-cbind(DM.CpGIs.stat,'cytoband'=cytobands.of.DM.cpgis)

