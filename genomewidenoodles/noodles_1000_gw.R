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

RangedData.from.df.simple<-function(df){
	RangedData(
		space=df$space, 
		ranges=IRanges
		(
			start=as.numeric(df$start),
			end=as.numeric(df$end)
		)
	)
}

noodle.length<-1000

noodles.1000.with.methylation.loaded<-FALSE
# we can the whole thing to noodles.1000.with.methylation.Rda
if(file.exists('noodles.1000.with.methylation.Rda'))
	if ('noodles.1000.with.methylation' %in% load('noodles.1000.with.methylation.Rda'))
		if (class(noodles.1000.with.methylation)=='data.frame')
			noodles.1000.with.methylation.loaded<-TRUE

if (!noodles.1000.with.methylation.loaded)
{
	source('../common/read_clinical.R')
	#Clinical prepared.

	#it is folder with bed files
	peakbedsfolder<-paste0(dataFolder,'/PeakCalls/bedfiles/')

	beds<-list.files(peakbedsfolder)

	#preparing noodles.1000
	# we can the whole thing to noodles.1000.with.methylation.Rda
	source('../common/load_or_read_chrom_ranges.R')

	noodles.space=character(0)
	noodles.start=integer(0)
	noodles.end=integer(0)
	chrom.length=end(ranges(chrom.ranges))
	message('generating noodles')
	for (chr_no in 1:length(chrom.ranges))
	{
		this.space<-paste('chr',as.character(seqnames(chrom.ranges)[chr_no]),sep='')
		message(this.space)
		names(chrom.length)[chr_no]<-this.space
		this.noodles.start<-seq(1,end(ranges(chrom.ranges))[chr_no],by=noodle.length)
		noodles.space<-c(noodles.space,rep(this.space,length(this.noodles.start)))
		noodles.start<-c(noodles.start,this.noodles.start)
		this.noodles.end<-this.noodles.start+noodle.length-1
		this.noodles.end[length(this.noodles.end)]<-min(this.noodles.end[length(this.noodles.end)],chrom.length[chr_no])
		noodles.end<-c(noodles.end,this.noodles.end)
	}
	message('combining noodles')
	noodles.1000.with.methylation<-data.frame(chr=noodles.space,start=noodles.start,end=noodles.end)
	noodles.1000<-makeGRangesFromDataFrame(noodles.1000.with.methylation,seqinfo=chrom.length)	

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
			methylcoverage.this.chr<-sapply(1:length(noodles.1000[chr][[1]]),function(band){
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
	message('Saving...')
	save(file='noodles.1000.with.methylation.Rda',list=c('noodles.1000.with.methylation','Clinical','clinFile','beds','bed_available','bed_used','tumors','normals','DNAids','noodle.length'))
}

#tests.number<-dim(noodles.1000.with.methylation)[1]

#wilcoxon.p.values<-numeric(tests.number)
#normals.are.less.methylated<-logical(tests.number)

#expected.w.statistic<-(sum(normals[bed_available])*sum(tumors[bed_available]))/2

#for (rown in 1:tests.number)
#{
#	meth.values<-jitter(as.numeric(noodles.1000.with.methylation[rown,][DNAids[bed_available]]))
#	w<-wilcox.test(meth.values[normals[bed_available]],meth.values[tumors[bed_available]])
#	wilcoxon.p.values[rown]<-w$p.value
#	normals.are.less.methylated[rown]<-(w[['statistic']]<expected.w.statistic)
#}

#wilcoxon.p.values.bonferroni<-p.adjust(wilcoxon.p.values,'bonferroni')
#wilcoxon.p.values.fdr<-p.adjust(wilcoxon.p.values,'fdr')

#DM.cytobands<-which(wilcoxon.p.values<=0.05)
#DM.cytobands.bonferroni<-which(wilcoxon.p.values.bonferroni<=0.05)
#DM.cytobands.fdr<-which(wilcoxon.p.values.fdr<=0.05)

