if (!suppressWarnings(require('rtracklayer')))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite("rtracklayer")
}

if (!suppressWarnings(require('DASiR')))
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
	source('../common/read_clinical.R')
	#Clinical prepared.

	#it is folder with bed files
	peakbedsfolder<-paste0(dataFolder,'/PeakCalls/bedfiles/')

	beds<-list.files(peakbedsfolder)

	#reading karyotype
	# we can the whole thing to karyotype.with.methylation.Rda
	source('../common/load_or_read_karyotype.R')
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
		message(DNAidKey)
		if (length(match)>1) stop(paste0("More than one match of DNAid ",DNAid," amonng the bed file names.\n"));
		bedfilename<-paste0(peakbedsfolder,beds[match[1]]);
		methylated.ranges<-as(import(bedfilename),"RangedData")
		overlaps<-findOverlaps(karyotype,methylated.ranges)
		methylcoverage=integer(0)
		for(chr in names(karyotype))
		#cycle by chromosome
		{
			list.of.ovelaps.in.this.chr<-as.list(overlaps[[chr]])
			width.of.meth.ranges.in.this.chr<-width(methylated.ranges[chr])
			methylcoverage.this.chr<-sapply(1:length(karyotype[chr][[1]]),function(band){
				#attention: length(CpGIs[chr][[1]] supposes that there is column in datarange other than space and ranges (e.g., Id)
				#otherwise, use  length(start(CpGIs[chr]))
				sum(width.of.meth.ranges.in.this.chr[list.of.ovelaps.in.this.chr[[band]]])
			})#list of methylated coverage per cytoband
			methylcoverage<-c(methylcoverage,methylcoverage.this.chr)
			#add to main list
			#we need dual cycle because of the findOverlap return structure
		}
		message('done\n')
		karyotype.with.methylation[[DNAid]]=methylcoverage
		bed_used[match[1]]<-TRUE
		bed_available<-c(bed_available,TRUE)
	}
	message('Saving...')
	save(file='karyotype.with.methylation.Rda',list=c('karyotype.with.methylation','Clinical','clinFile','beds','bed_available','bed_used','tumors','normals','DNAids'))
}

karyotype.DM.loaded<-FALSE
# we can the whole thing to karyotype.with.methylation.Rda
if(file.exists('karyotype.DM.Rda'))
	if ('DM.cytobands' %in% load('karyotype.DM.Rda'))
			karyotype.DM.loaded-TRUE

{
	tests.number<-dim(karyotype.with.methylation)[1]

	wilcoxon.p.values<-numeric(tests.number)
	normals.are.less.methylated<-logical(tests.number)

	expected.w.statistic<-(sum(normals[bed_available])*sum(tumors[bed_available]))/2

	for (rown in 1:tests.number)
	{
		meth.values<-jitter(as.numeric(karyotype.with.methylation[rown,][DNAids[bed_available]]))
		w<-wilcox.test(meth.values[normals[bed_available]],meth.values[tumors[bed_available]])
		wilcoxon.p.values[rown]<-w$p.value
		normals.are.less.methylated[rown]<-(w[['statistic']]<expected.w.statistic)
		#anova.result<-c(anova.result,anova(lm(meth.values~tumors))[1,'Pr(>F)'])
	}

	wilcoxon.p.values.bonferroni<-p.adjust(wilcoxon.p.values,'bonferroni')
	wilcoxon.p.values.fdr<-p.adjust(wilcoxon.p.values,'fdr')

	DM.cytobands<-which(wilcoxon.p.values<=0.05)
	DM.cytobands.bonferroni<-which(wilcoxon.p.values.bonferroni<=0.05)
	DM.cytobands.fdr<-which(wilcoxon.p.values.fdr<=0.05)
	save(file='karyotype.DM.Rda',list=c('DM.cytobands','DM.cytobands.bonferroni','DM.cytobands.fdr'))
}
