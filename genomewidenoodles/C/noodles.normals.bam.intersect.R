#we suppose all the preparatio is already done
library('Matrix')
library('Rsamtools')
library('GenomicAlignments')
if (!'noodles.C' %in% ls()) load('noodles.C.Rda')

BAM.folder<-'f:/Califano/MBDseq_PeakCalls/bam-bai_files/'

bamsinfolder<-list.files(BAM.folder,pattern='bam$')

normal.ids<-DNAids[contrast==0]

noodles.C.normals.read.coverage.loaded<-FALSE
# we can the whole thing to noodles.M.Rda
if(file.exists('noodles.C.normals.read.coverage.Rda'))
{
	loaded<-load('noodles.C.normals.read.coverage.Rda')
	if ('noodles.C.normals.read.coverage' %in% loaded) 
		if (class(noodles.C.normals.read.coverage)=='dgCMatrix' || 
				class(noodles.C.normals.read.coverage)=='matrix')
			noodles.C.normals.read.coverage.loaded<-TRUE
}

if(!noodles.C.normals.read.coverage.loaded)
{
	message('scanning bam file names')
	normal.bam.names<-unlist(lapply(normal.ids,function(id){
		id<-strsplit(id,',')[[1]][1]	#remove all after ,	
		sample_files<-bamsinfolder[grep(id,bamsinfolder)]
		enriched_sample_files<-sample_files[grep('nrich',sample_files)]
		}))
	message('done')

	noodles.C.normals.read.coverage<-Matrix(0,ncol=length(normal.ids),nrow=length(noodles.C),sparse = TRUE)

	colnames(noodles.C.normals.read.coverage)<-normal.ids

	for (sample.no in (1:length(normal.ids)))
	{
		message(normal.ids[sample.no])
		noodles.C.normals.read.coverage[,sample.no]<-countOverlaps(
			noodles.C,
			readGAlignments(
				paste0(BAM.folder,normal.bam.names[sample.no])
		))
	}
	save(file='noodles.C.normals.read.coverage.Rda',list=c('noodles.C.normals.read.coverage','normal.bam.names','BAM.folder'))
}
