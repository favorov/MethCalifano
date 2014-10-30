if (!suppressWarnings(require('differential.coverage')))
{
	if (!suppressWarnings(require('devtools')))
	{
		source("http://bioconductor.org/biocLite.R")
		biocLite("devtools")
		library("devtools")
	}
	install_github('favorov/Differential.Coverage')
	#load_all('../../../../differential.coverage/')
	library('differential.coverage')
}

noodles.M.loaded<-FALSE
# we can the whole thing to noodles.M.Rda
if(file.exists('noodles.M.Rda'))
{
	loaded<-load('noodles.M.Rda')
	if ('noodles.M.methylation' %in% loaded) 
		if (class(noodles.M.methylation)=='data.frame')
			if ('noodles.M' %in% loaded)
				if(class(noodles.M)=='GRanges')
			noodles.M.loaded<-TRUE
}

if(!noodles.M.loaded)
{
	beddir<-'../../../../Methylation/bedfiles/'
	noodle.length<-1000
	chrs<-nucl.chromosomes.hg19()
	noodles.M<-prepare.uniform.noodles(chrs,noodle.length)
	bedfiles<-dir(beddir) 
	bedfiles<-bedfiles[grep('All_',bedfiles,invert=TRUE)] # remove two 'All_' files
	bed.ids<-sapply(strsplit(bedfiles,split='_'),function(x){if(x[2]!='DNA') x[2] else x[3]}) #somhere id in pos 2, somewhere in 3
	noodles.M.methylation<-CountCoverageOfNoodles(noodles.M,paste0(beddir,bedfiles),bed.ids)
	save(file='noodles.M.Rda',list=c('noodles.M','noodles.M.methylation','bed.ids','noodle.length'))
}
