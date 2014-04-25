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

noodle.length<-100

noodles.100.with.wilcoxon.loaded<-FALSE
# we can the whole thing to noodles.100.with.wilcoxon.Rda
if(file.exists('noodles.100.with.wilcoxon.Rda'))
	if ('noodles.100.with.wilcoxon' %in% load('noodles.100.with.wilcoxon.Rda'))
		if (class(noodles.100.with.wilcoxon)=='data.frame')
			noodles.100.with.wilcoxon.loaded<-TRUE

if (!noodles.100.with.wilcoxon.loaded)
{
	
	noodles.100.with.methylation.loaded<-FALSE
	#the file name is noodles.100.with.methylation.Rda
	#it contains daframe noodles.100.df and methylation matrices
	#noodle.100.meth.in.normals and noodle.100.meth.in.tumors
	#as well as some work stuff. If we have it, we are not to 
	#rerun the noodler and the methylation interval overlapper
	if(file.exists('noodles.100.with.methylation.Rda'))
	{
		loaded<-load('noodles.100.with.methylation.Rda')
		if ( 
			('noodles.100.df' %in% loaded)  &&	(class(noodles.100.df)=='data.frame')
			&&
			('noodles.100.meth.in.normals'%in% loaded)  &&	(class(noodles.100.meth.in.normals)=='matrix')
			&&
			('noodles.100.meth.in.tumors'%in% loaded)  &&	(class(noodles.100.meth.in.normals)=='matrix')
			&&
			dim(noodles.100.df)[1]==dim(noodles.100.meth.in.tumors)[1]
			&&
			dim(noodles.100.df)[1]==dim(noodles.100.meth.in.normals)[1]
		)		noodles.100.with.methylation.loaded<-TRUE
	}
	
	if(!noodles.100.with.methylation.loaded)
	{
		source('../common/read_clinical.R')
		#Clinical prepared.

		#it is folder with bed files
		peakbedsfolder<-paste0(dataFolder,'/PeakCalls/bedfiles/')

		beds<-list.files(peakbedsfolder)

		#preparing noodles.100
		# we can the whole thing to noodles.100.with.methylation.Rda
		source('../common/load_or_read_chrom_ranges.R')

		noodles.100.space=character(0)
		noodles.100.start=integer(0)
		noodles.100.end=integer(0)
		chrom.length=end(ranges(chrom.ranges))
		message('generating noodles')
		for (chr_no in 1:length(chrom.ranges))
		{
			this.space<-paste('chr',as.character(seqnames(chrom.ranges)[chr_no]),sep='')
			message(this.space)
			names(chrom.length)[chr_no]<-this.space
			this.noodles.start<-seq(1,end(ranges(chrom.ranges))[chr_no],by=noodle.length)
			noodles.100.space<-c(noodles.100.space,rep(this.space,length(this.noodles.start)))
			noodles.100.start<-c(noodles.100.start,this.noodles.start)
			this.noodles.end<-this.noodles.start+noodle.length-1
			this.noodles.end[length(this.noodles.end)]<-min(this.noodles.end[length(this.noodles.end)],chrom.length[chr_no])
			noodles.100.end<-c(noodles.100.end,this.noodles.end)
		}
		message('combining noodles')
		noodles.100.df<-data.frame(chr=noodles.100.space,start=noodles.100.start,end=noodles.100.end)
		noodles.100.with.methylation<-noodles.100.df # we start with them same

		noodles.100<-RangedData(
			space=noodles.100.space, 
			ranges=IRanges
			(
				start=noodles.100.start,
				end=noodles.100.end
			)
		)

		message('noodles done')

		bed_available<-logical(0)
		bed_used<-rep(FALSE,length(beds))

		message('coverage..')

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
			overlaps<-findOverlaps(noodles.100,methylated.ranges)
			methylcoverage=integer(0)
			for(chr in names(noodles.100))
			#cycle by chromosome
			{
				list.of.ovelaps.in.this.chr<-as.list(overlaps[[chr]])
				width.of.meth.ranges.in.this.chr<-width(methylated.ranges[chr])
				methylcoverage.this.chr<-sapply(1:length(start(noodles.100[chr])),function(band){
					sum(width.of.meth.ranges.in.this.chr[list.of.ovelaps.in.this.chr[[band]]])
				})#list of methylated coverage per cytoband
				methylcoverage<-c(methylcoverage,methylcoverage.this.chr)
				#add to main list
				#we need dual cycle because of the findOverlap return structure
			}
			message('done\n')
			noodles.100.with.methylation[[DNAid]]=methylcoverage
			bed_used[match[1]]<-TRUE
			bed_available<-c(bed_available,TRUE)
		}

		message('Wilcoxon prep')
		noodles.100.meth<-noodles.100.with.methylation[,DNAids[bed_available]]
		noodles.100.meth.in.normals<-as.matrix(noodles.100.meth[,normals[bed_available]])
		noodles.100.meth.in.tumors<-as.matrix(noodles.100.meth[,tumors[bed_available]])

		message('done\n')
		message('Saving meth data')
		save(file='noodles.100.with.methylation.Rda',list=c('noodles.100.df','noodles.100.meth.in.normals','noodles.100.meth.in.tumors','Clinical','clinFile','beds','bed_available','bed_used','tumors','normals','DNAids','noodle.length'))
		message('done\n')
	
	}
	
	message('Wilcoxoning...')
	tests.number<-dim(noodles.100.df)[1]
	noodles.100.wilcoxon.p.values<-numeric(tests.number)
	noodles.100.normals.are.less.methylated<-logical(tests.number)

	expected.w.statistic<-(sum(normals[bed_available])*sum(tumors[bed_available]))/2

	for (rown in 1:tests.number)
	{
		if ((rown %% 10000)==0 ){message(paste(as.character(rown),'of',as.character(tests.number)))}
		if (max(noodles.100.meth.in.normals[rown,],noodles.100.meth.in.tumors[rown,])==0)
		{
				noodles.100.wilcoxon.p.values[rown]<-1
				next
		}
		#meth.values<-as.numeric(noodles.100.with.methylation[rown,][DNAids[bed_available]])
		#meth.values<-jitter(meth.values)
		#wilcoxor.res<-wilcox.test(meth.values[normals[bed_available]],meth.values[tumors[bed_available]])
		w<-wilcox.test(jitter(noodles.100.meth.in.normals[rown,]),jitter(noodles.100.meth.in.tumors[rown,]))
		noodles.100.wilcoxon.p.values[rown]<-w$p.value
		noodles.100.normals.are.less.methylated[rown]<-(w[['statistic']]<expected.w.statistic)
	}
	noodles.100.with.wilcoxon<-cbind(noodles.100.df,'p-value'=noodles.100.wilcoxon.p.values,'if.hyper'=noodles.100.normals.are.less.methylated)
	message('Saving...')
	save(file='noodles.100.with.wilcoxon.Rda',list=c('noodles.100.with.wilcoxon','Clinical','clinFile','beds','bed_available','bed_used','tumors','normals','DNAids','noodle.length'))
	message('done\n')
}

noodles.100.wilcoxon.p.values<-noodles.100.with.wilcoxon$'p-value'
noodles.100.wilcoxon.p.values.bonferroni<-p.adjust(noodles.100.wilcoxon.p.values,'bonferroni')
noodles.100.wilcoxon.p.values.fdr<-p.adjust(noodles.100.wilcoxon.p.values,'fdr')

DM.noodles.100.ids<-which(noodles.100.wilcoxon.p.values<=0.05)
DM.noodles.100.bonferroni.ids<-which(noodles.100.wilcoxon.p.values.bonferroni<=0.05)
DM.noodles.100.fdr.ids<-which(noodles.100.wilcoxon.p.values.fdr<=0.05)

DM.noodles.100.bonferroni<-noodles.100.with.wilcoxon[DM.noodles.100.bonferroni.ids,]




