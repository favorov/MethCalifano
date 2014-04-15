#combines uscs gened and HGNC ids. Also, provide promoters, 35 and 53 genes
#writes tables for direc and reverse genes  and for upstreams 
#saves the information to ../common/CpGIs.Rda
#if the file exists, read it after some checks


refseqGenesWithHGNC.loaded<-FALSE

if(file.exists('../common/refseqGenesWithHGNCids.rda'))
	if ('refseqWithHGNC' %in% load('../common/refseqGenesWithHGNCids.rda'))
		if (class(refseqWithHGNC)=='RangedData')
			refseqGenesWithHGNCids.loaded<-TRUE

if(!refseqGenesWithHGNC.loaded)
{

	source('../common/load_or_read_chrom_ranges.R')
	source('../common/load_or_read_HGNC_ids.R')
	source('../common/load_or_read_refseq_promoters.R')
	refseqGenesdf<-as(refseqGenes,'data.frame')
	promoterdf<-as(refseqPromoters,'data.frame')
	refseqTSSdf<-as(refseqTSS,'data.frame')

	refseq.hgnc.ids=hgnc.ids[grep('N[MR]',hgnc.ids$RefSeq.IDs),]
	#we save only those HGNC that have NM or NR notation
	#in Refseq, we have only NM and MR, so it's OK
	refseq.hgnc.ids$RefSeq.IDs[which(refseq.hgnc.ids$RefSeq.IDs=='NM_014686')[2]] = 'NM_014686_1'
	#it was the only duplication in MR,NM refseq ids in hgnc
	rownames(refseq.hgnc.ids)=refseq.hgnc.ids$RefSeq.IDs
	if.refseq.in.hgnc<-as.character(refseqGenesdf$label) %in% refseq.hgnc.ids$RefSeq.IDs
	refseqGenesWithHGNCidsdf<-refseqGenesdf[if.refseq.in.hgnc,]
	refseqGenesWithHGNCidsdf<-cbind(refseqGenesWithHGNCidsdf,refseq.hgnc.ids[as.character(refseqGenesWithHGNCidsdf$label),setdiff(colnames(refseq.hgnc.ids),'RefSeq.IDs')])
	#setdiff not to save RefSeq.IDs two times; actually, this invocation is a JOIN+SELECT: we addree refseq.hgnc.ids by the RefSeq.Ids and we join two tables


	if(0)
	{
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

	save(file='../common/refseqWithHGNC.rda',list=c('refseqPromoters','refseqGenes.35','refseqGenes.53','refseqTSS','promoter.after.start','promoter.before.start'))
	}
}


