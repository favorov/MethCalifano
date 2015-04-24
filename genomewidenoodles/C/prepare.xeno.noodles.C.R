#here, we just want to prepare 0/1 text table (noodle_number lines)
#for 4 samples: 
#2 xenos:
#1 is 43684 or 42132
#2 is 42133

#2 cell line models:
#SCC-090
#SCC-047-2

#then, perl script will converge the superhuge table with this table

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

xeno.C.methylation.loaded<-FALSE
# we can the whole thing to noodles.M.Rda
if(file.exists('xeno.C.methylation.Rda'))
{
	loaded<-load('xeno.C.methylation.Rda')
	if ('xeno.C.methylation' %in% loaded) 
		if (class(xeno.C.methylation)=='dgCMatrix' || 
				class(xeno.C.methylation)=='matrix')
			xeno.C.methylation.loaded<-TRUE
}

if(!xeno.C.methylation.loaded)
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
	source('../../common/prepare_beds_and_contrast.R')
	#beds and contrast is ready
	#it is folder with bed files

	xeno.C.methylation<-count.coverage.of.noodles(noodles.C,xbedfilenames,xDNAids)
	save(file='xeno.C.methylation.Rda',list=c('xeno.C.methylation'))
	
	#
	#this call prepares the file with matrix
}

