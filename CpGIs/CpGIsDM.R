if (!suppressWarnings(require('differential.coverage')))
{
	if (!suppressWarnings(require('devtools')))
	{
		source("http://bioconductor.org/biocLite.R")
		biocLite("devtools")
		library("devtools")
	}
	install_github('favorov/differential.coverage')
	#load_all('../../../../differential.coverage/')
	library('differential.coverage')
}

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
		if (class(CpGIs.methylation)=='dgCMatrix' || 
				class(CpGIs.methylation)=='matrix')
			CpGIs.methylation.loaded<-TRUE

if (!CpGIs.methylation.loaded)
{
	source('../common/read_clinical.R')
	#Clinical prepared.
	source('../common/prepare_beds_and_contrast.R')
	#beds and contrast is ready
	CpGIs.methylation<-count.coverage.of.noodles(CpGIs,bedfilenames,DNAids)
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

	set.seed(1248312)
	#to avoid different jittr result in different runs
	
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

#we want to put each diffmet CpGi to a cytoband

load('../CytoBands/cytobands.DM.Rda')


generate.DM.CpGi.report<-function(DM.CpGIslands.set,#indices
												set.id) #variable part of the output file names
{
	message('Generating report for ',set.id,'\n')
	tsvfilename=paste0("DM.CpGIs.stat.",set.id,".tsv")
	htmlfilename=paste0("DM.CpGIs.stat.",set.id,".html")

	report.file.info<-file.info(c(tsvfilename,htmlfilename))

	if ( !any(is.na(report.file.info$size)) && !any(report.file.info$size==0)) 
	{
		message(paste0('Both reports for ',set.id,' were already present; doing nothing\n'))
		return(NA)
	}
	
	DM.CpGIs.stat<-data.frame(
		'id'=elementMetadata(CpGIs)$id[DM.CpGIslands.set],
		'chr'=as.character(seqnames(CpGIs))[DM.CpGIslands.set],
		'start'=start(CpGIs)[DM.CpGIslands.set],
		'end'=end(CpGIs)[DM.CpGIslands.set],
		'wilcoxon.p.value'=wilcoxon.p.values[DM.CpGIslands.set],
		'hyper?'=normals.are.less.methylated[DM.CpGIslands.set],
		'fisher.p.value'=fisher.results$fisher.p.values[DM.CpGIslands.set],
		'tmr.ratio'=fisher.results$meth.in.tumors.ratio[DM.CpGIslands.set],
		'nor.ratio'=fisher.results$meth.in.normals.ratio[DM.CpGIslands.set],
		'OR'=fisher.results$OR[DM.CpGIslands.set],
		'CI_95_L'=fisher.results$CI_95_L[DM.CpGIslands.set],
		'CI_95_H'=fisher.results$CI_95_H[DM.CpGIslands.set]
	)

	rownames(DM.CpGIs.stat)<-NULL

	message('Mapping to karyotype...')

	DM.CpGIs<-CpGIs[DM.CpGIslands.set]

	CpGIs.to.karyotype<-findOverlaps(DM.CpGIs,cytobands,type="within")

	DM.CpGIs.cytobands<-sapply(1:length(DM.CpGIs),function(i)
		{
			cb.index<-which(i==queryHits(CpGIs.to.karyotype))
			#it could be 1 or 0, because of 'whthin'
			if (length(cb.index)==0) return (c(NA,NA))	
			cb<-subjectHits(CpGIs.to.karyotype)[cb.index]
			c(cb,(cytobands.DM.statistics$'wilcoxon.p.values'[cb]<0.05))
		}
	)

	DM.CpGIs.stat<-cbind(DM.CpGIs.stat,'cytoband'=cytobands$'name'[DM.CpGIs.cytobands[1,]],'DM.band?'=as.logical(DM.CpGIs.cytobands[2,]))
	message('done\n')

	message('Looking for closest genes')
	DM.CpGIs.closest.genes<-closest.gene.start.by.interval(DM.CpGIs)

	DM.CpGIs.stat<-cbind(DM.CpGIs.stat,elementMetadata(DM.CpGIs.closest.genes)[,c('closest.TSS','pos','dir','dist')])

	message('done')

	message('Looking for overlapped genes')

	flanks<-7000

	DM.CpGIs.ovelapped.genes<-genes.with.TSS.covered.by.interval(DM.CpGIs,flanks=flanks)

	DM.CpGIs.stat<-cbind(DM.CpGIs.stat,elementMetadata(DM.CpGIs.ovelapped.genes)[,c('overlapped.TSS','overlapped.pos','ovrl.dir')])

	message('done\n')

	DM.CpGIs.stat$id<-substr(DM.CpGIs.stat$id,6,1000) # 1000 'any'; we strip first 'CpGi: ' from the id

	#now, we order it according to GRanges order

	DM.CpGIs.stat<-DM.CpGIs.stat[order(DM.CpGIs),]
	
	write.table(DM.CpGIs.stat,file=tsvfilename,sep='\t',row.names=FALSE,quote=FALSE)

	if(file.exists(htmlfilename)) {file.remove(htmlfilename)}

	print(xtable(DM.CpGIs.stat,digits=c(0,0,0,0,0,8,0,8,2,2,2,2,2,0,0,0,0,0,0,0,0,0), display=c('d','s','s','d','d','g','s','g','f','f','f','f','f','s','s','s','d','s','d','s','s','s')), type="html", file=htmlfilename,include.rownames=FALSE)
	
	0
}

CpGIs.DM.indices.loaded<-FALSE
# we can the whole thing to CpGIs.methylation.Rda
if(file.exists('CpGIs.DM.indices.Rda'))
	if ('DM.CpGIslands' %in% load('CpGIs.DM.indices.Rda')) CpGIs.DM.indices.loaded<-TRUE

if(!CpGIs.DM.indices.loaded)
{
	#bonferroni 
	DM.W.CpGIslands.Bonferroni<-which(p.adjust(wilcoxon.p.values,method='bonferroni')<=0.05)
	DM.F.CpGIslands.Bonferroni<-which(p.adjust(fisher.results$fisher.p.values,method='bonferroni')<=0.05)
	DM.CpGIslands.Bonferroni<-sort(union(DM.W.CpGIslands.Bonferroni,DM.F.CpGIslands.Bonferroni))


	#fdr 
	DM.W.CpGIslands.FDR<-which(p.adjust(wilcoxon.p.values,method='fdr')<=0.1)
	DM.F.CpGIslands.FDR<-which(p.adjust(fisher.results$fisher.p.values,method='fdr')<=0.1)
	DM.CpGIslands.FDR<-sort(union(DM.W.CpGIslands.FDR,DM.F.CpGIslands.FDR))

	#uncorr
	DM.W.CpGIslands<-which(wilcoxon.p.values<=0.05)
	DM.F.CpGIslands<-which(fisher.results$fisher.p.values<=0.05)
	DM.CpGIslands<-sort(union(DM.W.CpGIslands,DM.F.CpGIslands))

	save(file='CpGIs.DM.indices.Rda',list=c('DM.CpGIslands.Bonferroni','DM.CpGIslands.FDR','DM.CpGIslands'))
}
#we generate the reports
message('Generating reports')
generate.DM.CpGi.report(DM.CpGIslands.Bonferroni,'bonf')
generate.DM.CpGi.report(DM.CpGIslands.FDR,'fdr')
generate.DM.CpGi.report(DM.CpGIslands,'uncorr')
message('done\n')
