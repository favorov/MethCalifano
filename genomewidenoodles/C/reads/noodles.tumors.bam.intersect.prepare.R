#we suppose all the preparatio is already done
BAM.folder<-'~/Califano/MBDseq_PeakCalls/bam-bai_files/'
#BAM.folder<-'g:/Califano/MBDseq_PeakCalls/bam-bai_files/'
noodles.Rda.folder<-'./'
noodles<-'noodles.C.7.spaghetti' 
noodles.bed.file<-paste0(noodles,'.bed')
# we suppose Rda is $noodles.Rda
noodles.Rda.file<-paste0(noodles.Rda.folder,noodles,'.Rda')
# the finale result is $noodles.tumorss.read.coverage
result<-paste0(noodles,'.tumors.read.coverage')
resultfilename<-paste0(result,'.Rda')

library('Matrix')
library('rtracklayer')

load('../noodles.C.Rda') #for dnaids, etc

if (!noodles %in% ls()) load (noodles.Rda.file)

bamsinfolder<-list.files(BAM.folder,pattern='bam$')

tumor.ids<-DNAids[contrast==1]

result.loaded<-FALSE
# we can the whole thing to noodles.M.Rda
if(file.exists(resultfilename))
{
	loaded<-load(resultfilename)
	if (result %in% loaded) 
		if (class(get(result))=='dgCMatrix' || 
				class(get(result))=='matrix')
			result.loaded<-TRUE
} 

if(!result.loaded)
{
	if (!(file.exists(noodles.bed.file) && (file.info(noodles.bed.file)$size>0)))
		export(get(noodles),con=noodles.bed.file,format='bed')

	message('scanning bam file names')
	tumor.bam.names<-unlist(lapply(tumor.ids,function(id){
		id<-strsplit(id,',')[[1]][1]	#remove all after ,	
		sample_files<-bamsinfolder[grep(id,bamsinfolder)]
		enriched_sample_files<-sample_files[grep('nrich',sample_files)]
		}))
	message('done')

	for (sample.no in 1:length(tumor.ids))
	{
		message(tumor.ids[sample.no])
		bamfilename<-tumor.bam.names[sample.no]
		bamfilefullpath<-paste0(BAM.folder,bamfilename)
		bambedfilename<-paste0(bamfilename,'.bed')
		countfilename<-paste0(bamfilename,'.count')
		bed.create.command<-paste0('bedtools bamtobed -i ',bamfilefullpath,' > ',bambedfilename)
		message(bed.create.command)
		intersect.command<-paste0('bedtools intersect -c -a ',noodles.bed.file,' -b ', bambedfilename, ' | cut -f7 > ',countfilename)
		scriptname<-paste0('intersect',tumor.ids[sample.no],'.sh')
		scriptname=gsub(', ','_',scriptname)
		sink(scriptname)
		cat('#!/bin/bash\n',bed.create.command,'\n',intersect.command,'\n',sep='')
		sink()
	}
}
