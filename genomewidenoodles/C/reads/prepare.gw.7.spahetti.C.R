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
#we form RDA file name for any spahetti as spaghetti name.Rda
noodle.C.7.spaghetti.loaded<-FALSE
# we can the whole thing to noodles.M.Rda
if(file.exists('noodles.C.7.spaghetti.Rda'))
{
	loaded<-load('noodles.C.7.spaghetti.Rda')
		if ('noodles.C.7.spaghetti' %in% loaded)
				if(class(noodles.C)=='GRanges')
			noodles.C.7.spaghetti.loaded<-TRUE
}

if(!noodles.C.7.spaghetti.loaded)
{
	noodle.length<-700
	chrs<-nucl.chromosomes.hg19()
	noodles.C.7.spaghetti<-prepare.covering.noodles(chrs,noodle.length)
	add.spagthetti<-GRanges()
	
	
	save(file='noodles.C.7.spaghetti.Rda',list=c('noodles.C.7.spaghetti.Rda','noodle.length'))
	
}

