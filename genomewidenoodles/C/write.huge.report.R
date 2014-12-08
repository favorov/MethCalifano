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

if(!('all.the.all.loaded' %in% ls()) || is.na(all.the.all.loaded)) all.the.all.loaded<-FALSE
#for quick-develop
if(!all.the.all.loaded)
{
	message('loading..')
	load('noodles.C.Rda')
	load('noodles.C.fisher.results.Rda')
	load('../../CytoBands/cytobands.DM.Rda')
	load('../../CpGIs/CpGIs.Rda')
	load('../../CpGIs/CpGIs.DM.indices.Rda')
	all.the.all.loaded<-TRUE
}

rows.no<-dim(fisher.noodles.C.result)[1]
report.interval<-1:rows.no

report.noodles<-noodles.C[report.interval,]
report.fisher<-fisher.noodles.C.result[report.interval,]
#actually, it is to develop for little tests.no

tsvfilename="noodles.C.complete.annotaion.tsv"
YesNofilename="noodles.C.methylation.tsv"


#prepare dataframe
message('init dataframe')
report.frame<-data.frame('chr'=as.character(seqnames(report.noodles)),start=start(report.noodles),end=end(report.noodles),stringsAsFactors = FALSE)

rownames(report.frame)<-paste(report.frame$chr,":",report.frame$start,'-',report.frame$end,sep='')

message('adding Fisher')
report.frame<-cbind(report.frame,
		'fisher.p.value'=report.fisher$fisher.p.values,
		'tmr.ratio'=report.fisher$meth.in.tumors.ratio,
		'nor.ratio'=report.fisher$meth.in.normals.ratio,
		'OR'=report.fisher$OR,
		'CI_95_L'=report.fisher$CI_95_L,
		'CI_95_H'=report.fisher$CI_95_H
	)

message('Mapping to karyotype...')

cb<-integer(rows.no)

noodles.to.karyotype<-findOverlaps(report.noodles,cytobands,type="within")

cb[queryHits(noodles.to.karyotype)]=subjectHits(noodles.to.karyotype)

cb[cb==0]=NA

message('done')

report.frame<-cbind(report.frame,'cytoband'=cytobands$'name'[cb],'DM.band?'=cytobands.DM.statistics$'wilcoxon.p.values'[cb]<0.05,stringsAsFactors = FALSE)
#prepared

message('Mapping to cpg islands...')

noodles.to.cpgi<-findOverlaps(report.noodles,CpGIs,type="within")

ci<-integer(rows.no)

ci[queryHits(noodles.to.cpgi)]=subjectHits(noodles.to.cpgi)

ci[ci==0]=NA

message('done')

report.frame<-cbind(report.frame,'CpGi'=CpGIs$'id'[ci],'DM.island?'=ifelse(is.na(ci),NA,as.logical(ci %in% DM.CpGIslands)),stringsAsFactors = FALSE)

#report.frame$'CpGi'<-substr(report.frame$'CpGi',6,1000) # 1000 'any'; we strip first 'CpGi: ' from the id

message('Looking for closest genes')
closest.genes<-closest.gene.start.by.interval(report.noodles)

message('done')

#message('Looking for overlapped genes')

#too loongloongloong time

#flanks<-7000

#ovelapped.genes<-genes.with.TSS.covered.by.interval(report.noodles,flanks=flanks)
#message('done')

message('combining')
report.frame<-cbind(report.frame,elementMetadata(closest.genes)[,c('closest.TSS','pos','dir','dist')],elementMetadata(ovelapped.genes)[,c('overlapped.TSS','overlapped.pos','ovrl.dir')])
report.frame<-cbind(report.frame,elementMetadata(closest.genes)[,c('closest.TSS','pos','dir','dist')])
message('done\n')

#prepared

message('Prepare the huge 0-1 matrix')
meth.mat<-ifelse(as.matrix(noodles.C.methylation[report.interval,]),1,0)
rownames(meth.mat)<-rownames(report.frame)
message('done')

save(file='huge.Rda',list=c('report.frame','meth.mat'))

message('merging')
report.frame<-cbind(report.frame,meth.mat)
message('done')

save(file='huge-merged.Rda',list=c('report.frame'))
#sink('noodles.C.complete.annotaion.tsv')

#sink()

#write.table(report.frame,file=tsvfilename,sep='\t',row.names=TRUE)

#write.table(met.mat,file=YesNofilename,sep='\t',row.names=TRUE,col.names=TRUE)

