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

source('../common/load_or_read_super_enhancers.R')


superEnhancers.methylation.loaded<-FALSE
# we can the whole thing to superEnhancers.methylation.Rda
if(file.exists('superEnhancers.methylation.Rda'))
	if ('superEnhancers.methylation' %in% load('superEnhancers.methylation.Rda'))
		if (class(superEnhancers.methylation)=='dgCMatrix' || 
				class(superEnhancers.methylation)=='matrix')
			superEnhancers.methylation.loaded<-TRUE


if (!superEnhancers.methylation.loaded)
{
	source('../common/read_clinical.R')
	#Clinical prepared.
	source('../common/prepare_beds_and_contrast.R')
	#beds and contrast is ready
	superEnhancers.methylation<-count.coverage.of.noodles(superEnhancers,bedfilenames,DNAids)
	#it is folder with bed files
	message('Saving...\n')
	save(file='superEnhancers.methylation.Rda',list=c('superEnhancers.methylation','Clinical','clinFile','clinFileName','bedsinfolder','bed.used','tumors','normals','contrast','DNAids'))
}

superEnhancers.wilcoxon.loaded<-FALSE
# we can the whole thing to superEnhancers.with.methylation.Rda
if(file.exists('superEnhancers.wilcoxon.Rda'))
	if ('wilcoxon.p.values' %in% load('superEnhancers.wilcoxon.Rda'))
			superEnhancers.wilcoxon.loaded<-TRUE

if(!superEnhancers.wilcoxon.loaded)
{
	message('Wilcoxon\n')

	tests.number<-dim(superEnhancers.methylation)[1]

	expected.w.statistic<-(sum(normals)*sum(tumors))/2

	wilcoxon.res<-apply(superEnhancers.methylation,1,function(row){
			meth.values<-jitter(row)
			w<-wilcox.test(meth.values[normals],meth.values[tumors])
			c(w$p.value,(w[['statistic']]<expected.w.statistic))
		})

	wilcoxon.p.values<-wilcoxon.res[1,]
	normals.are.less.methylated<-as.logical(wilcoxon.res[2,])

	message('Saving...\n')
	save(file='superEnhancers.wilcoxon.Rda',list=c('wilcoxon.p.values','normals.are.less.methylated','tests.number'))
	message('done\n')
}

superEnhancers.fisher.loaded<-FALSE
# we can the whole thing to superEnhancers.methylation.Rda
if(file.exists('superEnhancers.fisher.Rda'))
	if ('fisher.results' %in% load('superEnhancers.fisher.Rda') && 'data.frame' == class(fisher.results))
			superEnhancers.fisher.loaded<-TRUE


if(!superEnhancers.fisher.loaded)
{
	message('Fishering\n')
	tests.number<-dim(superEnhancers.methylation)[1]

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
		metraw<-superEnhancers.methylation[rown,]
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
	save(file='superEnhancers.fisher.Rda',list=c('fisher.results','tests.number','contrast'))
}

load('../CytoBands/cytobands.DM.Rda')

generate.DM.superEnhaners.report<-function(DM.superEnhancers.set,#indices
												set.id) #variable part of the output file names
{
	message('Generating report for ',set.id,'\n')
	
	DM.superEnhancers.stat<-data.frame(
		'chr'=as.character(seqnames(superEnhancers))[DM.superEnhancers.set],
		'start'=start(superEnhancers)[DM.superEnhancers.set],
		'end'=end(superEnhancers)[DM.superEnhancers.set],
		'shortlist'=elementMetadata(superEnhancers)$short.list[DM.superEnhancers.set],
		'type'=elementMetadata(superEnhancers)$type[DM.superEnhancers.set],
		'wilcoxon.p.value'=wilcoxon.p.values[DM.superEnhancers.set],
		'hyper?'=normals.are.less.methylated[DM.superEnhancers.set],
		'fisher.p.value'=fisher.results$fisher.p.values[DM.superEnhancers.set],
		'tmr.ratio'=fisher.results$meth.in.tumors.ratio[DM.superEnhancers.set],
		'nor.ratio'=fisher.results$meth.in.normals.ratio[DM.superEnhancers.set],
		'OR'=fisher.results$OR[DM.superEnhancers.set],
		'CI_95_L'=fisher.results$CI_95_L[DM.superEnhancers.set],
		'CI_95_H'=fisher.results$CI_95_H[DM.superEnhancers.set]
	)

	tsvfilename=paste0("DM.superEnhancers.stat.",set.id,".tsv")
	htmlfilename=paste0("DM.superEnhancers.stat.",set.id,".html")

	report.file.info<-file.info(c(tsvfilename,htmlfilename))

	if ( !any(is.na(report.file.info$size)) && !any(report.file.info$size==0)) 
	{
		message(paste0('Both reports for ',set.id,' were already present; doing nothing\n'))
		return(NA)
	}
	
	rownames(DM.superEnhancers.stat)<-NULL

	message('Mapping to karyotype...')
	
	DM.superEnhancers<-superEnhancers[DM.superEnhancers.set]

	superEnhancers.to.karyotype<-findOverlaps(DM.superEnhancers,cytobands,type="within")

	DM.superEnhancers.cytobands<-
		sapply(1:length(DM.superEnhancers),function(i)
			{
				cb.index<-which(i==queryHits(superEnhancers.to.karyotype))
				#it could be 1 or 0, because of 'whthin'
				if (length(cb.index)==0) return (c(NA,NA))	
				cb<-subjectHits(superEnhancers.to.karyotype)[cb.index]
				c(cb,(cytobands.DM.statistics$'wilcoxon.p.values'[cb]<0.05))
			}
		)

	DM.superEnhancers.stat<-cbind(DM.superEnhancers.stat,'cytoband'=cytobands$'name'[DM.superEnhancers.cytobands[1,]],'DM.band?'=as.logical(DM.superEnhancers.cytobands[2,]))
	message('done\n')

	message('Looking for closest genes')
	DM.superEnhancers.closest.genes<-closest.gene.start.by.interval(DM.superEnhancers)

	DM.superEnhancers.stat<-cbind(DM.superEnhancers.stat,elementMetadata(DM.superEnhancers.closest.genes)[,c('closest.TSS','pos','dir','dist')])

	message('done')

	message('Looking for overlapped genes')

	flanks<-7000

	DM.superEnhancers.ovelapped.genes<-genes.with.TSS.covered.by.interval(DM.superEnhancers,flanks=flanks)

	DM.superEnhancers.stat<-cbind(DM.superEnhancers.stat,elementMetadata(DM.superEnhancers.ovelapped.genes)[,c('overlapped.TSS','overlapped.pos','ovrl.dir')])

	message('done\n')

	#now, we order it according to GRanges order

	DM.superEnhancers.stat<-DM.superEnhancers.stat[order(DM.superEnhancers),]
	

	write.table(DM.superEnhancers.stat,file=tsvfilename,sep='\t',row.names=FALSE,quote=FALSE)


	if(file.exists(htmlfilename)) {file.remove(htmlfilename)}

	print(xtable(DM.superEnhancers.stat,digits=c(0,0,0,0,0,0,8,0,8,2,2,2,2,2,0,0,0,0,0,0,0,0,0), display=c('s','s','d','d','s','s','g','s','g','f','f','f','f','f','s','s','s','d','s','d','s','s','s')), type="html", file=htmlfilename,include.rownames=FALSE)

	0
}


superEnhancers.DM.indices.loaded<-FALSE
# we can the whole thing to superEnhancers.methylation.Rda
if(file.exists('superEnhancers.DM.indices.Rda'))
	if ('DM.superEnhancers' %in% load('superEnhancers.DM.indices.Rda')) superEnhancers.DM.indices.loaded<-TRUE

if(!superEnhancers.DM.indices.loaded)
{
	#bonferroni 
	DM.W.superEnhancers.Bonferroni<-which(p.adjust(wilcoxon.p.values,method='bonferroni')<=0.05)
	DM.F.superEnhancers.Bonferroni<-which(p.adjust(fisher.results$fisher.p.values,method='bonferroni')<=0.05)
	DM.superEnhancers.Bonferroni<-sort(union(DM.W.superEnhancers.Bonferroni,DM.F.superEnhancers.Bonferroni))
	DM.superEnhancers.Bonferroni.and<-sort(intersect(DM.W.superEnhancers.Bonferroni,DM.F.superEnhancers.Bonferroni))


	#fdr 
	DM.W.superEnhancers.FDR<-which(p.adjust(wilcoxon.p.values,method='fdr')<=0.1)
	DM.F.superEnhancers.FDR<-which(p.adjust(fisher.results$fisher.p.values,method='fdr')<=0.1)
	DM.superEnhancers.FDR<-sort(union(DM.W.superEnhancers.FDR,DM.F.superEnhancers.FDR))
	DM.superEnhancers.FDR.and<-sort(intersect(DM.W.superEnhancers.FDR,DM.F.superEnhancers.FDR))

	#uncorr
	DM.W.superEnhancers<-which(wilcoxon.p.values<=0.05)
	DM.F.superEnhancers<-which(fisher.results$fisher.p.values<=0.05)
	DM.superEnhancers<-sort(union(DM.W.superEnhancers,DM.F.superEnhancers))
	DM.superEnhancers.and<-sort(intersect(DM.W.superEnhancers,DM.F.superEnhancers))

	save(file='superEnhancers.DM.indices.Rda',
		list=c(
			'DM.F.superEnhancers.Bonferroni',
			'DM.F.superEnhancers.FDR',
			'DM.F.superEnhancers',
			'DM.W.superEnhancers.Bonferroni',
			'DM.W.superEnhancers.FDR',
			'DM.W.superEnhancers',
			'DM.superEnhancers.Bonferroni',
			'DM.superEnhancers.FDR',
			'DM.superEnhancers',
			'DM.superEnhancers.Bonferroni.and',
			'DM.superEnhancers.FDR.and',
			'DM.superEnhancers.and'))
}
#we generate the reports
message('Generating reports')
generate.DM.superEnhaners.report(DM.superEnhancers.Bonferroni,'bonf')
generate.DM.superEnhaners.report(DM.superEnhancers.FDR,'fdr')
generate.DM.superEnhaners.report(DM.superEnhancers,'uncorr')
generate.DM.superEnhaners.report(DM.superEnhancers.Bonferroni.and,'bonf.and')
generate.DM.superEnhaners.report(DM.superEnhancers.FDR.and,'fdr.and')
generate.DM.superEnhaners.report(DM.superEnhancers.and,'uncorr.and')
message('done\n')

