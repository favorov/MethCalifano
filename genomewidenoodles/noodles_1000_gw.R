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

noodle.length<-1000

noodles.1000.with.wilcoxon.loaded<-FALSE
# we can the whole thing to noodles.1000.with.wilcoxon.Rda
if(file.exists('noodles.1000.with.wilcoxon.Rda'))
	if ('noodles.1000.with.wilcoxon' %in% load('noodles.1000.with.wilcoxon.Rda'))
		if (class(noodles.1000.with.wilcoxon)=='data.frame')
			noodles.1000.with.wilcoxon.loaded<-TRUE

if (!noodles.1000.with.wilcoxon.loaded)
{
	
	noodles.1000.with.methylation.loaded<-FALSE
	#the file name is noodles.1000.with.methylation.Rda
	#it contains daframe noodles.1000.df and methylation matrices
	#noodle.1000.meth.in.normals and noodle.1000.meth.in.tumors
	#as well as some work stuff. If we have it, we are not to 
	#rerun the noodler and the methylation interval overlapper
	if(file.exists('noodles.1000.with.methylation.Rda'))
	{
		loaded<-load('noodles.1000.with.methylation.Rda')
		if ( 
			('noodles.1000.df' %in% loaded)  &&	(class(noodles.1000.df)=='data.frame')
			&&
			('noodles.1000.meth.in.normals'%in% loaded)  &&	(class(noodles.1000.meth.in.normals)=='matrix')
			&&
			('noodles.1000.meth.in.tumors'%in% loaded)  &&	(class(noodles.1000.meth.in.normals)=='matrix')
			&&
			dim(noodles.1000.df)[1]==dim(noodles.1000.meth.in.tumors)[1]
			&&
			dim(noodles.1000.df)[1]==dim(noodles.1000.meth.in.normals)[1]
		)		noodles.1000.with.methylation.loaded<-TRUE
	}
	
	if(!noodles.1000.with.methylation.loaded)
	{
		source('../common/read_clinical.R')
		#Clinical prepared.

		#it is folder with bed files
		peakbedsfolder<-paste0(dataFolder,'/PeakCalls/bedfiles/')

		beds<-list.files(peakbedsfolder)

		#preparing noodles.1000
		# we can the whole thing to noodles.1000.with.methylation.Rda
		source('../common/load_or_read_chrom_ranges.R')

		noodles.1000.space=character(0)
		noodles.1000.start=integer(0)
		noodles.1000.end=integer(0)
		chrom.length=end(ranges(chrom.ranges))
		message('generating noodles')
		for (chr_no in 1:length(chrom.ranges))
		{
			this.space<-paste('chr',as.character(seqnames(chrom.ranges)[chr_no]),sep='')
			message(this.space)
			names(chrom.length)[chr_no]<-this.space
			this.noodles.start<-seq(1,end(ranges(chrom.ranges))[chr_no],by=noodle.length)
			noodles.1000.space<-c(noodles.1000.space,rep(this.space,length(this.noodles.start)))
			noodles.1000.start<-c(noodles.1000.start,this.noodles.start)
			this.noodles.end<-this.noodles.start+noodle.length-1
			this.noodles.end[length(this.noodles.end)]<-min(this.noodles.end[length(this.noodles.end)],chrom.length[chr_no])
			noodles.1000.end<-c(noodles.1000.end,this.noodles.end)
		}
		message('combining noodles')
		noodles.1000.df<-data.frame(chr=noodles.1000.space,start=noodles.1000.start,end=noodles.1000.end)
		noodles.1000.with.methylation<-noodles.1000.df # we start with them same

		noodles.1000<-RangedData(
			space=noodles.1000.space, 
			ranges=IRanges
			(
				start=noodles.1000.start,
				end=noodles.1000.end
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
			overlaps<-findOverlaps(noodles.1000,methylated.ranges)
			methylcoverage=integer(0)
			for(chr in names(noodles.1000))
			#cycle by chromosome
			{
				list.of.ovelaps.in.this.chr<-as.list(overlaps[[chr]])
				width.of.meth.ranges.in.this.chr<-width(methylated.ranges[chr])
				methylcoverage.this.chr<-sapply(1:length(start(noodles.1000[chr])),function(band){
					sum(width.of.meth.ranges.in.this.chr[list.of.ovelaps.in.this.chr[[band]]])
				})#list of methylated coverage per cytoband
				methylcoverage<-c(methylcoverage,methylcoverage.this.chr)
				#add to main list
				#we need dual cycle because of the findOverlap return structure
			}
			message('done\n')
			noodles.1000.with.methylation[[DNAid]]=methylcoverage
			bed_used[match[1]]<-TRUE
			bed_available<-c(bed_available,TRUE)
		}

		message('Wilcoxon prep')
		noodles.1000.meth<-noodles.1000.with.methylation[,DNAids[bed_available]]
		noodles.1000.meth.in.normals<-as.matrix(noodles.1000.meth[,normals[bed_available]])
		noodles.1000.meth.in.tumors<-as.matrix(noodles.1000.meth[,tumors[bed_available]])

		message('done\n')
		message('Saving meth data')
		save(file='noodles.1000.with.methylation.Rda',list=c('noodles.1000.df','noodles.1000.meth.in.normals','noodles.1000.meth.in.tumors','Clinical','clinFile','beds','bed_available','bed_used','tumors','normals','DNAids','noodle.length'))
		message('done\n')
	
	}
	
	message('Wilcoxoning...')
	tests.number<-dim(noodles.1000.df)[1]
	noodles.1000.wilcoxon.p.values<-numeric(tests.number)
	noodles.1000.normals.are.less.methylated<-logical(tests.number)

	expected.w.statistic<-(sum(normals[bed_available])*sum(tumors[bed_available]))/2

	for (rown in 1:tests.number)
	{
		if ((rown %% 100000)==0 ){message(paste(as.character(rown),'of',as.character(tests.number)))}
		if (max(noodles.1000.meth.in.normals[rown,],noodles.1000.meth.in.tumors[rown,])==0)
		{
				noodles.1000.wilcoxon.p.values[rown]<-1
				next
		}
		#meth.values<-as.numeric(noodles.1000.with.methylation[rown,][DNAids[bed_available]])
		#meth.values<-jitter(meth.values)
		#wilcoxor.res<-wilcox.test(meth.values[normals[bed_available]],meth.values[tumors[bed_available]])
		w<-wilcox.test(jitter(noodles.1000.meth.in.normals[rown,]),jitter(noodles.1000.meth.in.tumors[rown,]))
		noodles.1000.wilcoxon.p.values[rown]<-w$p.value
		noodles.1000.normals.are.less.methylated[rown]<-(w[['statistic']]<expected.w.statistic)
	}
	noodles.1000.with.wilcoxon<-cbind(noodles.1000.df,'p-value'=noodles.1000.wilcoxon.p.values,'if.hyper'=noodles.1000.normals.are.less.methylated)
	message('Saving...')
	save(file='noodles.1000.with.wilcoxon.Rda',list=c('noodles.1000.with.wilcoxon','Clinical','clinFile','beds','bed_available','bed_used','tumors','normals','DNAids','noodle.length'))
	message('done\n')
}

noodles.1000.wilcoxon.p.values<-noodles.1000.with.wilcoxon$'p-value'
noodles.1000.wilcoxon.p.values.bonferroni<-p.adjust(noodles.1000.wilcoxon.p.values,'bonferroni')
noodles.1000.wilcoxon.p.values.fdr<-p.adjust(noodles.1000.wilcoxon.p.values,'fdr')

DM.noodles.1000.ids<-which(noodles.1000.wilcoxon.p.values<=0.05)
DM.noodles.1000.bonferroni.ids<-which(noodles.1000.wilcoxon.p.values.bonferroni<=0.05)
DM.noodles.1000.fdr.ids<-which(noodles.1000.wilcoxon.p.values.fdr<=0.05)

DM.noodles.1000.bonferroni<-noodles.1000.with.wilcoxon[DM.noodles.1000.bonferroni.ids,]




