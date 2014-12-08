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

source('../common/load_or_read_vista_enhancers.R')


vistaEnhancers.methylation.loaded<-FALSE
# we can the whole thing to vistaEnhancers.methylation.Rda
if(file.exists('vistaEnhancers.methylation.Rda'))
	if ('vistaEnhancers.methylation' %in% load('vistaEnhancers.methylation.Rda'))
		if (class(vistaEnhancers.methylation)=='dgCMatrix' || 
				class(vistaEnhancers.methylation)=='matrix')
			vistaEnhancers.methylation.loaded<-TRUE


if (!vistaEnhancers.methylation.loaded)
{
	source('../common/read_clinical.R')
	#Clinical prepared.
	source('../common/prepare_beds_and_contrast.R')
	#beds and contrast is ready
	vistaEnhancers.methylation<-count.coverage.of.noodles(vistaEnhancers,bedfilenames,DNAids)
	#it is folder with bed files
	message('Saving...\n')
	save(file='vistaEnhancers.methylation.Rda',list=c('vistaEnhancers.methylation','Clinical','clinFile','clinFileName','bedsinfolder','bed.used','tumors','normals','contrast','DNAids'))
}

vistaEnhancers.wilcoxon.loaded<-FALSE
# we can the whole thing to vistaEnhancers.with.methylation.Rda
if(file.exists('vistaEnhancers.wilcoxon.Rda'))
	if ('wilcoxon.p.values' %in% load('vistaEnhancers.wilcoxon.Rda'))
			vistaEnhancers.wilcoxon.loaded<-TRUE

if(!vistaEnhancers.wilcoxon.loaded)
{
	message('Wilcoxon\n')

	tests.number<-dim(vistaEnhancers.methylation)[1]

	expected.w.statistic<-(sum(normals)*sum(tumors))/2

	wilcoxon.res<-apply(vistaEnhancers.methylation,1,function(row){
			meth.values<-jitter(row)
			w<-wilcox.test(meth.values[normals],meth.values[tumors])
			c(w$p.value,(w[['statistic']]<expected.w.statistic))
		})

	wilcoxon.p.values<-wilcoxon.res[1,]
	normals.are.less.methylated<-as.logical(wilcoxon.res[2,])

	message('Saving...\n')
	save(file='vistaEnhancers.wilcoxon.Rda',list=c('wilcoxon.p.values','normals.are.less.methylated','tests.number'))
	message('done\n')
}

vistaEnhancers.fisher.loaded<-FALSE
# we can the whole thing to vistaEnhancers.methylation.Rda
if(file.exists('vistaEnhancers.fisher.Rda'))
	if ('fisher.results' %in% load('vistaEnhancers.fisher.Rda') && 'data.frame' == class(fisher.results))
			vistaEnhancers.fisher.loaded<-TRUE


if(!vistaEnhancers.fisher.loaded)
{
	message('Fishering\n')
	tests.number<-dim(vistaEnhancers.methylation)[1]

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
		metraw<-vistaEnhancers.methylation[rown,]
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
	save(file='vistaEnhancers.fisher.Rda',list=c('fisher.results','tests.number','contrast'))
}

load('../CytoBands/cytobands.DM.Rda')

generate.DM.vistaEnhaners.report<-function(DM.vistaEnhancers.set,#indices
												set.id) #variable part of the output file names
{
	message('Generating report for ',set.id,'\n')
	
	DM.vistaEnhancers.stat<-data.frame(
		'id'=elementMetadata(vistaEnhancers)$id[DM.vistaEnhancers.set],
		'chr'=as.character(seqnames(vistaEnhancers))[DM.vistaEnhancers.set],
		'start'=start(vistaEnhancers)[DM.vistaEnhancers.set],
		'end'=end(vistaEnhancers)[DM.vistaEnhancers.set],
		'wilcoxon.p.value'=wilcoxon.p.values[DM.vistaEnhancers.set],
		'hyper?'=normals.are.less.methylated[DM.vistaEnhancers.set],
		'fisher.p.value'=fisher.results$fisher.p.values[DM.vistaEnhancers.set],
		'tmr.ratio'=fisher.results$meth.in.tumors.ratio[DM.vistaEnhancers.set],
		'nor.ratio'=fisher.results$meth.in.normals.ratio[DM.vistaEnhancers.set],
		'OR'=fisher.results$OR[DM.vistaEnhancers.set],
		'CI_95_L'=fisher.results$CI_95_L[DM.vistaEnhancers.set],
		'CI_95_H'=fisher.results$CI_95_H[DM.vistaEnhancers.set]
	)

	tsvfilename=paste0("DM.vistaEnhancers.stat.",set.id,".tsv")
	htmlfilename=paste0("DM.vistaEnhancers.stat.",set.id,".html")

	report.file.info<-file.info(c(tsvfilename,htmlfilename))

	if ( !any(is.na(report.file.info$size)) && !any(report.file.info$size==0)) 
	{
		message(paste0('Both reports for ',set.id,' were already present; doing nothing\n'))
		return(NA)
	}
	
	rownames(DM.vistaEnhancers.stat)<-NULL

	message('Mapping to karyotype...')
	
	DM.vistaEnhancers<-vistaEnhancers[DM.vistaEnhancers.set]

	vistaEnhancers.to.karyotype<-findOverlaps(DM.vistaEnhancers,cytobands,type="within")

	DM.vistaEnhancers.cytobands<-sapply(1:length(DM.vistaEnhancers),function(i)
		{
			cb<-subjectHits(vistaEnhancers.to.karyotype)[which(i==queryHits(vistaEnhancers.to.karyotype))]
			c(cb,(cytobands.DM.statistics$'wilcoxon.p.values'[cb]<0.05))
		}
	)

	DM.vistaEnhancers.stat<-cbind(DM.vistaEnhancers.stat,'cytoband'=cytobands$'name'[DM.vistaEnhancers.cytobands[1,]],'DM.band?'=as.logical(DM.vistaEnhancers.cytobands[2,]))
	message('done\n')

	message('Looking for closest genes')
	DM.vistaEnhancers.closest.genes<-closest.gene.start.by.interval(DM.vistaEnhancers)

	DM.vistaEnhancers.stat<-cbind(DM.vistaEnhancers.stat,elementMetadata(DM.vistaEnhancers.closest.genes)[,c('closest.TSS','pos','dir','dist')])

	message('done')

	message('Looking for overlapped genes')

	flanks<-7000

	DM.vistaEnhancers.ovelapped.genes<-genes.with.TSS.covered.by.interval(DM.vistaEnhancers,flanks=flanks)

	DM.vistaEnhancers.stat<-cbind(DM.vistaEnhancers.stat,elementMetadata(DM.vistaEnhancers.ovelapped.genes)[,c('overlapped.TSS','overlapped.pos','ovrl.dir')])

	message('done\n')

	DM.vistaEnhancers.stat$id<-substr(DM.vistaEnhancers.stat$id,6,1000) # 1000 'any'; we strip first 'CpGi: ' from the id

	#now, we order it according to GRanges order

	DM.vistaEnhancers.stat<-DM.vistaEnhancers.stat[order(DM.vistaEnhancers),]
	

	write.table(DM.vistaEnhancers.stat,file=tsvfilename,sep='\t',row.names=FALSE)


	if(file.exists(htmlfilename)) {file.remove(htmlfilename)}

	print(xtable(DM.vistaEnhancers.stat,digits=c(0,0,0,0,0,8,0,8,2,2,2,2,2,0,0,0,0,0,0,0,0,0), display=c('d','s','s','d','d','g','s','g','f','f','f','f','f','s','s','s','d','s','d','s','s','s')), type="html", file=htmlfilename,include.rownames=FALSE)
	
	0
}


vistaEnhancers.DM.indices.loaded<-FALSE
# we can the whole thing to vistaEnhancers.methylation.Rda
if(file.exists('vistaEnhancers.DM.indices.Rda'))
	if ('DM.vistaEnhancers' %in% load('vistaEnhancers.DM.indices.Rda')) vistaEnhancers.DM.indices.loaded<-TRUE

if(!vistaEnhancers.DM.indices.loaded)
{
	#bonferroni 
	DM.W.vistaEnhancers.Bonferroni<-which(p.adjust(wilcoxon.p.values,method='bonferroni')<=0.05)
	DM.F.vistaEnhancers.Bonferroni<-which(p.adjust(fisher.results$fisher.p.values,method='bonferroni')<=0.05)
	DM.vistaEnhancers.Bonferroni<-sort(union(DM.W.vistaEnhancers.Bonferroni,DM.F.vistaEnhancers.Bonferroni))
	DM.vistaEnhancers.Bonferroni.and<-sort(intersect(DM.W.vistaEnhancers.Bonferroni,DM.F.vistaEnhancers.Bonferroni))


	#fdr 
	DM.W.vistaEnhancers.FDR<-which(p.adjust(wilcoxon.p.values,method='fdr')<=0.1)
	DM.F.vistaEnhancers.FDR<-which(p.adjust(fisher.results$fisher.p.values,method='fdr')<=0.1)
	DM.vistaEnhancers.FDR<-sort(union(DM.W.vistaEnhancers.FDR,DM.F.vistaEnhancers.FDR))
	DM.vistaEnhancers.FDR.and<-sort(intersect(DM.W.vistaEnhancers.FDR,DM.F.vistaEnhancers.FDR))

	#uncorr
	DM.W.vistaEnhancers<-which(wilcoxon.p.values<=0.05)
	DM.F.vistaEnhancers<-which(fisher.results$fisher.p.values<=0.05)
	DM.vistaEnhancers<-sort(union(DM.W.vistaEnhancers,DM.F.vistaEnhancers))
	DM.vistaEnhancers.and<-sort(intersect(DM.W.vistaEnhancers,DM.F.vistaEnhancers))

	save(file='vistaEnhancers.DM.indices.Rda',
		list=c(
			'DM.F.vistaEnhancers.Bonferroni',
			'DM.F.vistaEnhancers.FDR',
			'DM.F.vistaEnhancers',
			'DM.W.vistaEnhancers.Bonferroni',
			'DM.W.vistaEnhancers.FDR',
			'DM.W.vistaEnhancers',
			'DM.vistaEnhancers.Bonferroni',
			'DM.vistaEnhancers.FDR',
			'DM.vistaEnhancers',
			'DM.vistaEnhancers.Bonferroni.and',
			'DM.vistaEnhancers.FDR.and',
			'DM.vistaEnhancers.and'))
}
#we generate the reports
message('Generating reports')
generate.DM.vistaEnhaners.report(DM.vistaEnhancers.Bonferroni,'bonf')
generate.DM.vistaEnhaners.report(DM.vistaEnhancers.FDR,'fdr')
generate.DM.vistaEnhaners.report(DM.vistaEnhancers,'uncorr')
generate.DM.vistaEnhaners.report(DM.vistaEnhancers.Bonferroni.and,'bonf.and')
generate.DM.vistaEnhaners.report(DM.vistaEnhancers.FDR.and,'fdr.and')
generate.DM.vistaEnhaners.report(DM.vistaEnhancers.and,'uncorr.and')
message('done\n')

