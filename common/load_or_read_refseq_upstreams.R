#combines uscs gened and HGNC ids.
#writes tables for direc and reverse genes  and for upstreams 
#saves the information to ../common/CpGIs.Rda
#if the file exists, read it after some checks

source('../common/load_or_read_chrom_ranges.R')
source('../common/load_or_read_known_genes.R')
#source('../common/load_or_read_HGNC_ids.R')

promoter.before.start=1000
promoter.after.start=0

refseqGenesdf<-as(refseqGenes,'data.frame')
genes.35<-refseqGenesdf[['orientation']]=='-'
genes.53<-refseqGenesdf[['orientation']]=='+'
refseqGenesdf.35<-refseqGenesdf[genes.35,]
refseqGenesdf.53<-refseqGenesdf[genes.53,]
promoterdf<-data.frame(
	apply(refseqGenesdf,1,
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
		)
,stringsAsFactors=FALSE)

#for (row.n in 1:nrow(promoterdf))
#{
#	if(promoterdf[row.n,'orientation']=='+')
#	{
#		promoterdf[row.n,'end']=min(as.integer(promoterdf[row.n,'start'])+promoter.after.start,chrom.length[promoterdf[row.n,'space']])
#		promoterdf[row.n,'start']=max(as.integer(promoterdf[row.n,'start'])-promoter.before.start,1)
#	}
#	else
#	{
#		promoterdf[row.n,'start']=max(as.integer(promoterdf[row.n,'end'])-promoter.after.start,1)
#		promoterdf[row.n,'end']=min(as.integer(promoterdf[row.n,'end'])+promoter.before.start,chrom.length[promoterdf[row.n,'space']])
#	}
#	promoterdf[row.n,'width']=promoterdf[row.n,'end']-promoterdf[row.n,'start']+1
#}


