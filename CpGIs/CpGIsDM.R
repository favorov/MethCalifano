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

source('../common/load_or_read_CpGIs.R')

CpGIs.methylation.loaded<-FALSE
# we can the whole thing to CpGIs.methylation.Rda
if(file.exists('CpGIs.methylation.Rda'))
	if ('CpGIs.methylation' %in% load('CpGIs.methylation.Rda'))
		if (class(CpGIs.methylation)=='data.frame')
			CpGIs.methylation.loaded<-TRUE

if (!CpGIs.methylation.loaded)
{
	source('../common/read_clinical.R')
	#Clinical prepared.
	source('../common/prepare_beds_and_contrast.R ')
	#beds and contrast is ready
	CpGIs.methylation<-CountCoverageOfNoodles(CpGIs,bedfilenames,DNAids)
	#it is folder with bed files
	message('Saving...\n')
	save(file='CpGIs.methylation.Rda',list=c('CpGIs.methylation','Clinical','clinFile','clinFileName','bedsinfolder','bed.used','tumors','normals','contrast','DNAids'))
}

CpGIs.wilcoxon.loaded<-FALSE
# we can the whole thing to CpGIs.with.methylation.Rda
if(file.exists('CpGIs.wilcoxon.Rda'))
	if ('wilcoxon.p.values' %in% load('CpGIs.wilcoxon.Rda'))
			CpGIs.wilcoxon.loaded<-TRUE

if(!CpGIs.wilcoxon.loaded)
{
	message('Wilcoxon\n')

	tests.number<-dim(CpGIs.methylation)[1]

	expected.w.statistic<-(sum(normals)*sum(tumors))/2

	wilcoxon.res<-apply(CpGIs.methylation,1,function(row){
			meth.values<-jitter(row)
			w<-wilcox.test(meth.values[normals],meth.values[tumors])
			c(w$p.value,(w[['statistic']]<expected.w.statistic))
		})

	wilcoxon.p.values<-wilcoxon.res[1,]
	normals.are.less.methylated<-as.logical(wilcoxon.res[2,])

	message('Saving...\n')
	save(file='CpGIs.wilcoxon.Rda',list=c('wilcoxon.p.values','normals.are.less.methylated','tests.number'))
	message('done\n')
}


CpGIs.fisher.loaded<-FALSE
# we can the whole thing to CpGIs.methylation.Rda
if(file.exists('CpGIs.fisher.Rda'))
	if ('fisher.results' %in% load('CpGIs.fisher.Rda') && 'data.frame' == class(fisher.results))
			CpGIs.fisher.loaded<-TRUE


if(!CpGIs.fisher.loaded)
{
	message('Fishering\n')
	tests.number<-dim(CpGIs.methylation)[1]

	norm.no<-length(which(normals))
	tumor.no<-length(which(tumors))

	fishtabs<-as.matrix(prepare.tabulated.fisher(tumor.no,norm.no))

	message('create result matrix')
	fisher.noodles.result.mat<-matrix(fishtabs[1,],ncol=6,nrow=tests.number,byrow=TRUE)
	
	colnames(fisher.noodles.result.mat)<-c('fisher.p.values','meth.in.normals.ratio','meth.in.tumors.ratio','OR','CI_95_L','CI_95_H') 
	
	revcontrast<-!contrast
	report.every<-tests.number %/% 100
	message('fill result matrix')

	for (rown in 1:tests.number) 	
	{
		if (!(rown %% report.every)) message(rown)
		metraw<-CpGIs.methylation[rown,]
		aslogic<-as.logical(metraw)
		MY<-sum(aslogic & contrast)
		MN<-sum(aslogic & revcontrast)
		if (0==MN && 0==MY) next
		fishres<-fishtabs[tab.fisher.row.no(tumor.no,norm.no,MY,MN),]
		fisher.noodles.result.mat[rown,]<-fishres
	}

	message('converting to dataframe')
	fisher.results<-as.data.frame(fisher.noodles.result.mat)
	message('done\n')
	message('Saving...\n')
	save(file='CpGIs.fisher.Rda',list=c('fisher.results','tests.number','contrast'))
}

browser()
#we want to put each diffmet CpGi to a cytoband

load('../CytoBands/karyotype.methylation.Rda')
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

