if (!suppressWarnings(require('Matrix')))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite("Matrix")
	library("Matrix")
}
load('noodles.C.7.spaghetti.normals.read.coverage.Rda')

#if (!suppressWarnings(require('caTools')))
#{
#	source("http://bioconductor.org/biocLite.R")
#	biocLite("caTools")
#	library("caTools")
#}
#spaghetti.size.in.noodles<-7
#spaghetti.C.normals.read.coverage<-
#	spaghetti.size.in.noodles*
#	caTools::runmean(noodles.C.normals.read.coverage,spaghetti.size.in.noodles,alg='fast')
#running mean*window.size is running sum
#S4Vectors::runmean tries to shade the caTools::runmean

norm.read.stats.frame<-t(apply(noodles.C.7.spaghetti.normals.read.coverage,1,quantile))
save(file='noodles.C.7.spaghetti.normals.read.quantiles.Rda',list=c('norm.read.stats.frame'))
