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

if (!require('Homo.sapiens'))
{
  source("http://bioconductor.org/biocLite.R")
  biocLite("Homo.sapiens")
  library('Homo.sapiens')  
}
if (!require('org.Hs.eg.db'))
{
  source("http://bioconductor.org/biocLite.R")
  biocLite('org.Hs.eg.db')
  library('org.Hs.eg.db')  
}

if(!require('xtable'))
{
  source("http://bioconductor.org/biocLite.R")
  biocLite("xtable")
  library("xtable")  
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
	source('../common/load_or_read_CpGIs.R')
	#islands are read fro DAS or loaded

	CpGIs.with.methylation<-as(CpGIs,"data.frame")

	bed_available<-logical(0)
	bed_used<-rep(FALSE,length(beds))

	for (DNAid in DNAids)
	{ 
		DNAidKey<-strsplit(DNAid,',')[[1]][1]	#remove all after ,	
		message(DNAidKey)
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
				list.of.ovelaps.in.this.chr<-as.list(overlaps[[chr]])
				width.of.meth.ranges.in.this.chr<-width(methylated.ranges[chr])
				methylcoverage.this.chr<-sapply(1:length(CpGIs[chr][[1]]),function(island_no){
					#attention: length(CpGIs[chr][[1]] supposes that there is column in datarange other than space and ranges (e.g., Id)
					#otherwise, use  length(start(CpGIs[chr]))
					sum(width.of.meth.ranges.in.this.chr[list.of.ovelaps.in.this.chr[[island_no]]])
				})#list of methylated coverage per cytoband
			#list of methylated coverage per island
			#the main idea is that as.list(hitsobject) convert is to list of vectors, 
			#list is addressd by query oblects and give ref object lists (possibly, empty)
			methylcoverage<-c(methylcoverage,methylcoverage.this.chr)
			#add to main list
			#we need dual cycle because of the findOverlap return structure
		}
		CpGIs.with.methylation[[DNAid]]=methylcoverage
		bed_used[match[1]]<-TRUE
		bed_available<-c(bed_available,TRUE)
		message('done\n')
	}
	message('Saving...\n')
	save(file='CpGIs.with.methylation.Rda',list=c('CpGIs.with.methylation','Clinical','clinFile','beds','bed_available','bed_used','tumors','normals','DNAids'))
}

CpGIs.wilcoxon.data.loaded<-FALSE
# we can the whole thing to CpGIs.with.methylation.Rda
if(file.exists('CpGIs.wilcoxon.data.Rda'))
	if ('wilcoxon.p.values' %in% load('CpGIs.wilcoxon.data.Rda'))
			CpGIs.wilcoxon.data.loaded<-TRUE
if(!CpGIs.wilcoxon.data.loaded)
{
	message('Wilcoxon\n')
	wilcoxon.p.values<-numeric(0)

	normals.are.less.methyl.covered<-logical(0)

	expected.w.statistic<-(sum(normals[bed_available])*sum(tumors[bed_available]))/2

	tests.number<-dim(CpGIs.with.methylation)[1]

	for (rown in 1:tests.number)
	{
		meth.values<-jitter(as.numeric(CpGIs.with.methylation[rown,][DNAids[bed_available]]))
		w<-wilcox.test(meth.values[normals[bed_available]],meth.values[tumors[bed_available]])
		wilcoxon.p.values<-c(wilcoxon.p.values,w$p.value)
		normals.are.less.methyl.covered<-c(normals.are.less.methyl.covered,(w[['statistic']]<expected.w.statistic))
		#anova.result<-c(anova.result,anova(lm(meth.values~tumors))[1,'Pr(>F)'])
	}

	message('done\n')
	message('Saving...\n')
	save(file='CpGIs.wilcoxon.data.Rda',list=c('wilcoxon.p.values','normals.are.less.methyl.covered','tests.number'))
}

CpGIs.fisher.data.loaded<-FALSE
# we can the whole thing to CpGIs.with.methylation.Rda
if(file.exists('CpGIs.fisher.data.Rda'))
	if ('fisher.p.values' %in% load('CpGIs.fisher.data.Rda'))
			CpGIs.fisher.data.loaded<-TRUE
if(!CpGIs.fisher.data.loaded)
{
	message('Fishering\n')
	tests.number<-dim(CpGIs.with.methylation)[1]
	fisher.p.values<-numeric(tests.number)
	meth.in.normals.ratio<-numeric(tests.number)
	meth.in.tumors.ratio<-numeric(tests.number)
	OR<-numeric(tests.number)
	CI_95_L<-numeric(tests.number)
	CI_95_H<-numeric(tests.number)


	for (rown in 1:tests.number)
	{
		cotable<-table(as.logical(as.numeric(CpGIs.with.methylation[rown,][DNAids[bed_available]])),tumors[bed_available])
		if(nrow(cotable)==1)#nonmeth
		{
			fisher.p.values[rown]<-1.
			meth.in.tumors.ratio[rown]<-0
			meth.in.normals.ratio[rown]<-0
			OR[rown]<-NA
			CI_95_L[rown]<-NA
			CI_95_H[rown]<-NA
			next
		}
		fisherres<-fisher.test(cotable)
		fisher.p.values[rown]<-fisherres$p.value
		meth.in.tumors.ratio[rown]<-cotable[2,2]/cotable[1,2]
		meth.in.normals.ratio[rown]<-cotable[2,1]/cotable[1,1]
		OR[rown]<-fisherres$estimate
		CI_95_L[rown]<-fisherres$conf.int[1]
		CI_95_H[rown]<-fisherres$conf.int[2]
	}

	message('done\n')
	message('Saving...\n')
	save(file='CpGIs.fisher.data.Rda',list=c('fisher.p.values','tests.number','meth.in.tumors.ratio','meth.in.normals.ratio','OR','CI_95_L','CI_95_H'))
}

#we want to put each diffmet CpGi to a cytoband
source('../common/load_or_read_karyotype.R')
load('../CytoBands/karyotype.with.methylation.Rda')
load('../CytoBands/karyotype.DM.Rda')



generate.DM.CpGi.report<-function(DM.CpGIslands.set,#indices
												set.id) #variable part of the output file names
{
	message('Generating report for ',set.id,'\n')
	columns<-c('id','space','start','end')
	DM.CpGIs.stat<-cbind(
		CpGIs.with.methylation[DM.CpGIslands.set,columns],
		'wilcoxon.p.value'=wilcoxon.p.values[DM.CpGIslands.set],
		'hyper?'=normals.are.less.methyl.covered[DM.CpGIslands.set],
		'fisher.p.value'=fisher.p.values[DM.CpGIslands.set],
		'tmr.ratio'=meth.in.tumors.ratio[DM.CpGIslands.set],
		'nor.ratio'=meth.in.normals.ratio[DM.CpGIslands.set],
		'OR'=OR[DM.CpGIslands.set],
		'CI_95_L'=CI_95_L[DM.CpGIslands.set],
		'CI_95_H'=CI_95_H[DM.CpGIslands.set]
	)

	rownames(DM.CpGIs.stat)<-NULL

	message('Mapping to karyotype...')

	DM.CpGIs.Ranges<-as(DM.CpGIs.stat[,columns],'RangedData')


	CpGIs.to.karyotype<-findOverlaps(DM.CpGIs.Ranges,karyotype,type="within")
	cytobands.of.DM.cpgis=character(0)
	is.cytobands.of.DM.cpgis.DM=logical(0)

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
		cytobands.of.DM.cpgis.this.chr<-as.character(karyotype[chr][[1]][cytoband.numbers.this.chr])
		cytobands.of.DM.cpgis<-c(cytobands.of.DM.cpgis,cytobands.of.DM.cpgis.this.chr)
		is.cytobands.of.DM.cpgis.DM.this.chr<-cytobands.of.DM.cpgis.this.chr %in% as.character(karyotype.with.methylation$id)[intersect(which(karyotype.with.methylation$space==chr),DM.cytobands)]
		is.cytobands.of.DM.cpgis.DM<-c(is.cytobands.of.DM.cpgis.DM,is.cytobands.of.DM.cpgis.DM.this.chr)
	}


	DM.CpGIs.stat<-cbind(DM.CpGIs.stat,'cytoband'=cytobands.of.DM.cpgis,'DM.band?'=is.cytobands.of.DM.cpgis.DM)
	message('done\n')


	DM.CpGIs.GRanges<-GRanges(
		seqinfo=nucl.chromosomes.hg19(),
		ranges=IRanges(start=DM.CpGIs.stat[,'start'],end=DM.CpGIs.stat[,'end']),
		seqnames=DM.CpGIs.stat[,'space'],
		id=DM.CpGIs.stat[,'id']
	)

	seqinfo(DM.CpGIs.GRanges)<-nucl.chromosomes.hg19()

	message('Looking for closest genes')
	DM.CpGIs.closest.genes<-closest.gene.start.by.interval(DM.CpGIs.GRanges)

	DM.CpGIs.stat<-cbind(DM.CpGIs.stat,elementMetadata(DM.CpGIs.closest.genes)[,c('closest.TSS','pos','dir','dist')])

	message('done')

	message('Looking for overlapped genes')

	flanks<-7000

	DM.CpGIs.ovelapped.genes<-genes.with.TSS.covered.by.interval(DM.CpGIs.GRanges,flanks=flanks)

	DM.CpGIs.stat<-cbind(DM.CpGIs.stat,elementMetadata(DM.CpGIs.ovelapped.genes)[,c('overlapped.TSS','overlapped.pos','ovrl.dir')])

	message('done\n')

	DM.CpGIs.stat$id<-substr(DM.CpGIs.stat$id,6,1000) # 1000 'any'; we strip first 'CpGi: ' from the id

	hg.chr.order<-function(chrnames)
	{
		suppressWarnings(
			sapply(chrnames,function(chr)
				{
					if(!is.na(as.integer(as.character(chr)))) 
						return (as.integer(as.character(chr)))
					if(!is.na(as.integer(substring(as.character(chr),4)))) 
						return (as.integer(substring(as.character(chr),4)))
					if(chr=='X' || chr=='chrX') 
						return (22)
					if(chr=='Y' || chr=='chrY') 
						return (23)
					if(chr=='M' || chr=='chrM' || chr=='MT'|| chr=='chrMT') 
						return (24)
					return(NA)
				}
		))
	}

	DM.CpGIs.stat<-DM.CpGIs.stat[order(hg.chr.order(DM.CpGIs.stat$space),DM.CpGIs.stat$start),]
	
	tsvfilename=paste0("DM.CpGIs.stat.",set.id,".tsv")

	write.table(DM.CpGIs.stat,file=tsvfilename,sep='\t',row.names=FALSE)

	htmlfilename=paste0("DM.CpGIs.stat.",set.id,".html")

	if(file.exists(htmlfilename)) {file.remove(htmlfilename)}

	print(xtable(DM.CpGIs.stat,digits=c(0,0,0,0,0,8,0,8,2,2,2,2,2,0,0,0,0,0,0,0,0,0), display=c('d','s','s','d','d','g','s','g','f','f','f','f','f','s','s','s','d','s','d','s','s','s')), type="html", file=htmlfilename,include.rownames=FALSE)
	
	0
}

#we generate the reports
#bonferroni 
DM.W.CpGIslands.Bonferroni<-which(p.adjust(wilcoxon.p.values,method='bonferroni')<=0.05)
DM.F.CpGIslands.Bonferroni<-which(p.adjust(fisher.p.values,method='bonferroni')<=0.05)
DM.CpGIslands.Bonferroni<-sort(union(DM.W.CpGIslands.Bonferroni,DM.F.CpGIslands.Bonferroni))
generate.DM.CpGi.report(DM.CpGIslands.Bonferroni,'bonf')

#fdr 
DM.W.CpGIslands.FDR<-which(p.adjust(wilcoxon.p.values,method='fdr')<=0.1)
DM.F.CpGIslands.FDR<-which(p.adjust(fisher.p.values,method='fdr')<=0.1)
DM.CpGIslands.FDR<-sort(union(DM.W.CpGIslands.FDR,DM.F.CpGIslands.FDR))
generate.DM.CpGi.report(DM.CpGIslands.FDR,'fdr')

#uncorr
DM.W.CpGIslands<-which(wilcoxon.p.values<=0.05)
DM.F.CpGIslands<-which(fisher.p.values<=0.05)
DM.CpGIslands<-sort(union(DM.W.CpGIslands,DM.F.CpGIslands))
generate.DM.CpGi.report(DM.CpGIslands,'uncorr')

