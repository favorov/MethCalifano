#combines uscs gened and HGNC ids.
#writes tables for direc and reverse genes  and for upstreams 
#saves the information to ../common/CpGIs.Rda
#if the file exists, read it after some checks

source('../common/load_or_read_chrom_ranges.R')
source('../common/load_or_read_HGNC_ids.R')

refseqGenesWithHGNC.loaded<-FALSE

if(file.exists('../common/refseqGenesWithHGNCids.rda'))
	if ('refseqWithHGNC' %in% load('../common/refseqGenesWithHGNCids.rda'))
		if (class(refseqWithHGNC)=='RangedData')
			refseqPromoters.loaded<-TRUE

if(!refseqWithHGNC.loaded)
{

	refseqGenesdf<-as(refseqGenes,'data.frame')


	promoterdf<-data.frame(
		t(apply(refseqGenesdf,1,
				function(x)
				{
					if(x['orientation']=='+')
					{
						x['end']=min(as.integer(x['start'])+promoter.after.start,chrom.length[x['space']])
						x['start']=max(as.integer(x['start'])-promoter.before.start,1)
					}
					else
					{
						x['start']=max(as.integer(x['end'])-promoter.after.start,1)
						x['end']=min(as.integer(x['end'])+promoter.before.start,chrom.length[x['space']])
					}
					x['width']=NA
					x
				}
			))
	,stringsAsFactors=FALSE)

	refseqPromoters<-RangedData.from.df.refseq(promoterdf)

	refseqTSSdf<-data.frame(
		t(apply(refseqGenesdf,1,
				function(x)
				{
					if(x['orientation']=='+')
						x['end']=x['start']
					else
						x['start']=x['end']
					x['width']=NA
					x
				}
			))
	,stringsAsFactors=FALSE)

	refseqTSS<-RangedData.from.df.refseq(refseqTSSdf)

	save(file='../common/refseqPromoters.rda',list=c('refseqPromoters','refseqGenes.35','refseqGenes.53','refseqTSS','promoter.after.start','promoter.before.start'))
}


