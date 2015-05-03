#we suppose all the preparatio is already done
BAM.folder<-'~/Califano/MBDseq_PeakCalls/bam-bai_files/'
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
library('data.table')

message('loading noodles')
load('../noodles.C.Rda') #for dnaids, etc
message('done')

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
		id<-strsplit(id,'_')[[1]][1]	#remove all after _	
		sample_files<-bamsinfolder[grep(id,bamsinfolder)]
		enriched_sample_files<-sample_files[grep('nrich',sample_files)]
		}))
	message('done')

 	resultmatrix<-data.table(empty<-rep(0,length(noodles.C.7.spaghetti)))
	#Matrix(nrow=length(get(noodles)),ncol=0,sparse=TRUE)

	#resultmatrix<-matrix(0,ncol=length(tumor.ids),nrow=length(get(noodles)),sparse = TRUE)


	for (sample.no in 1:length(tumor.ids))
	{
		tumor.id<-tumor.ids[sample.no]
		message(tumor.ids[sample.no])
		bamfilename<-tumor.bam.names[sample.no]
		countfilename<-paste0(bamfilename,'.count')
		if(!(file.exists(countfilename) && file.info(countfilename)$size>0))
		{
			bamfilefullpath<-paste0(BAM.folder,bamfilename)
			bambedfilename<-paste0(bamfilename,'.bed')
			bed.create.command<-paste0('bedtools bamtobed -i ',bamfilefullpath,' > ',bambedfilename)
			message(bed.create.command)
			shell(bed.create.command)
			intersect.command<-paste0('bedtools intersect -c -a ',noodles.bed.file,' -b ', bambedfilename, ' | cut -f7 > ',countfilename)
			message(intersect.command)
			shell(intersect.command)
		}
		message('Reading...')
		data<-fread(countfilename)
		message('Binding...')
		resultmatrix[,eval(as.name(tumor.id)):=data]
		#unlink(c(bambedfilename,countfilename))
	}
	resultmatrix[,'V1':=NULL]
	assign(result,resultmatrix) # ref-copy
	save(file=resultfilename,list=c(result,'noodles','tumor.bam.names','BAM.folder'))
}
