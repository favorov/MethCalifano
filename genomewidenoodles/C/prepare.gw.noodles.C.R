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
	#it is folder with bed files
	peakbedsfolder<-paste0(data.folder,'/PeakCalls/bedfiles/')

	beds<-list.files(peakbedsfolder)

	DNAids<-DNAids[normals | tumors]
	#we do not want to work with xeongraft & cell lines

	contrast<-integer(0)

	contrast[normals & ! tumors] <- 0
	contrast[tumors & ! normals] <- 1
	#contrast : 0 for normals, 1 for tumors, NA for unknown 
	
	contrast<-contrast[!is.na(contrast)]
	#we do not want to work with xeongraft & cell lines

	bedfilenames<-lapply(DNAids,function(DNAid){ 
			DNAidKey<-strsplit(DNAid,',')[[1]][1]	#remove all after ,	
			message(DNAidKey)
			match<-grep(DNAidKey,beds)
			if (!length(match)) 
			{
				DNAidKey<-paste0(strsplit(DNAid,'_')[[1]],collapse='') 
				#remove _ from key; sometimes, it help
				match<-grep(DNAidKey,beds)
			}
			if (!length(match)) 
			{
				bed_available<-c(bed_available,FALSE)
				next
			}
			if (length(match)>1) stop(paste0("More than one match of DNAid ",DNAid," amonng the bed file names.\n"));
			bedfilename<-paste0(peakbedsfolder,beds[match[1]]);
			#message(match[1])
			#message(bedfilename)
			bedfilename
		}
	)
	message('Files assigned.\n')
	bedfilenames<-unlist(bedfilenames)
	noodles.C.methylation<-CountCoverageOfNoodles(noodles.C,bedfilenames,DNAids)
	save(file='noodles.C.Rda',list=c('noodles.C','noodles.C.methylation','DNAids','bedfilenames','contrast','noodle.length'))
	
	#
	norm.no<-length(which(contrast==0))
	tumor.no<-length(which(contrast==1))

	prepare.tabulated.fisher(tumor.no,norm.no)
	#this call prepares the file with matrix
}

