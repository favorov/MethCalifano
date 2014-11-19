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


chromosomes.meth.coverage.loaded<-FALSE
# we can the whole thing to cytobands.DM.Rda
if(file.exists('chromosomes.meth.coverage.Rda'))
	if ('chrom.meth.cover' %in% load('chromosomes.meth.coverage.Rda'))
		if (class(chrom.meth.cover)=='dgCMatrix' || 
				class(chrom.meth.cover)=='matrix')
			chromosomes.meth.coverage.loaded<-TRUE

chromosomes.meth.coverage.loaded<-FALSE

if (!chromosomes.meth.coverage.loaded)
{
	source('../common/read_clinical.R')
	#Clinical prepared.
	source('../common/prepare_beds_and_contrast.R')
	#bedfiles are ready
	chrinfo<-nucl.chromosomes.hg19()
	chrranges<-GRanges(seqinfo=chrinfo,ranges=IRanges(start=rep(1,24),end=seqlengths(chrinfo)),seqnames=seqnames(chrinfo))
	chrom.meth.cover<-as.matrix(CountCoverageOfNoodles(chrranges,bedfilenames,DNAids))
	rownames(chrom.meth.cover)<-seqnames(chrinfo)
	chrom.meth.cover.rel<-apply(chrom.meth.cover,2,function(cov) cov/seqlengths(chrinfo))
	save(file='chromosomes.meth.coverage.Rda',list=c('chrom.meth.cover','chrom.meth.cover.rel','chrranges','Clinical','clinFile','clinFileName','bedsinfolder','bed.used','tumors','normals','contrast','DNAids'))
}
