#combines uscs gened and HGNC ids.
#writes tables for direc and reverse genes  and for upstreams 
#saves the information to ../common/CpGIs.Rda
#if the file exists, read it after some checks

source('../common/load_or_read_chrom_ranges.R')
source('../common/load_or_read_refseq_genes.R')
#source('../common/load_or_read_HGNC_ids.R')

refseqPromoters.loaded<-FALSE

if(file.exists('../common/refseqPromoters.rda'))
	if ('refseqTSS' %in% load('../common/refseqPromoters.rda'))
		if (class(refseqTSS)=='RangedData')
			refseqPromoters.loaded<-TRUE

if(!refseqPromoters.loaded)
{
	promoter.before.start=1000
	promoter.after.start=0

	refseqGenesdf<-as(refseqGenes,'data.frame')

	genes.35<-refseqGenesdf[['orientation']]=='-'
	genes.53<-refseqGenesdf[['orientation']]=='+'

	refseqGenes.35<-RangedData.from.df.refseq(refseqGenesdf[genes.35,])
	refseqGenes.53<-RangedData.from.df.refseq(refseqGenesdf[genes.53,])

	refseqPromoterdf<-data.frame(
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

	refseqPromoters<-RangedData.from.df.refseq(refseqPromoterdf)

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


