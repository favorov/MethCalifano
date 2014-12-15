if(file.exists('noodles.C.7.spaghetti.Rda'))
{
	message('the file exists')
} else
{

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

	flanks<-300

	if (!'noodles.C' %in% ls()) load('../noodles.C.Rda')

	noodles.C.7.spaghetti<-noodles.C

	#inflate spahetti 
	start(noodles.C.7.spaghetti)<-pmax(1,start(noodles.C)-flanks)
	end(noodles.C.7.spaghetti)<-pmin(end(noodles.C)+flanks,as.integer(seqlengths(noodles.C)[as.character(seqnames(noodles.C))]))
	#inflated

	save(file='noodles.C.7.spaghetti.Rda',list=c('noodles.C.7.spaghetti'))
}		

