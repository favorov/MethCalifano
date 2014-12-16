#we suppose all the preparatio is already done
BAM.folder<-'~/Califano/MBDseq_PeakCalls/bam-bai_files/'
noodles.Rda.folder<-'./'
noodles<-'noodles.C.7.spaghetti' 
noodles.bed.file<-paste0(noodles,'.bed')
# we suppose Rda is $noodles.Rda
noodles.Rda.file<-paste0(noodles.Rda.folder,noodles,'.Rda')
# the finale result is $noodles.normals.read.coverage
result<-paste0(noodles,'.normals.read.coverage')
resultfilename<-paste0(result,'.Rda')

library('Matrix')
library('rtracklayer')

load('../noodles.C.Rda') #for dnaids, etc

if (!noodles %in% ls()) load (noodles.Rda.file)

bamsinfolder<-list.files(BAM.folder,pattern='bam$')

normal.ids<-DNAids[contrast==0]

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
		export(get(noodles,con=noodles.bed.file,format='bed'))

	message('scanning bam file names')
	normal.bam.names<-unlist(lapply(normal.ids,function(id){
		id<-strsplit(id,',')[[1]][1]	#remove all after ,	
		sample_files<-bamsinfolder[grep(id,bamsinfolder)]
		enriched_sample_files<-sample_files[grep('nrich',sample_files)]
		}))
	message('done')

	resultmatrix<-Matrix(0,ncol=length(normal.ids),nrow=length(get(noodles)),sparse = TRUE)

	colnames(resultmatrix)<-normal.ids

	for (sample.no in 1:length(normal.ids))
	{
		message(normal.ids[sample.no])
		bamfilename<-normal.bam.names[sample.no]
		bamfilefullpath<-paste0(BAM.folder,bamfilename)
		bambedfilename<-paste0(bamfilename,'.bed')
		countfilename<-paste0(bamfilename,'.count')
		bed.create.command<-paste0('bedtools bamtobed -i ',bamfilefullpath,' > ',bambedfilename)
		message(bed.create.command)
		shell(bed.create.command)
		intersect.command<-paste0('bedtools intersect -c -a ',noodles.bed.file,' -b ', bambedfilename, ' | cut -f7 > ',countfilename)
		message(intersect.command)
		shell(intersect.command)
		resultmatrix[,sample.no]=unlist(read.table(countfilename))
		unlink(c(bambedfilename,countfilename))
	}
	assign(result,resultmatrix)
	save(file=resultfilename,list=c(result,'noodles','normal.bam.names','BAM.folder'))
}
