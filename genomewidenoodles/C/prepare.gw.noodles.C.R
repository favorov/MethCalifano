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

noodles.C.loaded<-FALSE
# we can the whole thing to noodles.M.Rda
if(file.exists('noodles.C.Rda'))
{
	loaded<-load('noodles.C.Rda')
	if ('noodles.C.methylation' %in% loaded) 
		if (class(noodles.C.methylation)=='data.frame')
			if ('noodles.C' %in% loaded)
				if(class(noodles.C)=='GRanges')
			noodles.C.loaded<-TRUE
}

if(!noodles.C.loaded)
{
	noodle.length<-100
	chrs<-nucl.chromosomes.hg19()
	noodles.C<-prepare.covering.noodles(chrs,noodle.length)
	#noodles.C<-prepare.covering.noodles(chrs['chr1'],noodle.length)
	#chr1 test
	
	#project-dependent part
	#we read the IDS and codes from the clinical file
	data.version<-'final'
	data.folder<-'../../../Data'
	source('../../common/read_clinical.R')
	#Clinical prepared.
	source('../../common/prepare_beds_and_contrast.R ')
	#beds and contrast is ready
	#it is folder with bed files

	noodles.C.methylation<-CountCoverageOfNoodles(noodles.C,bedfilenames,DNAids)
	save(file='noodles.C.Rda',list=c('noodles.C','noodles.C.methylation','DNAids','bedfilenames','contrast','noodle.length'))
	
	#
	norm.no<-length(which(contrast==0))
	tumor.no<-length(which(contrast==1))

	prepare.tabulated.fisher(tumor.no,norm.no)
	#this call prepares the file with matrix
}

