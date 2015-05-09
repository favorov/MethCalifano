if (!suppressWarnings(require('data.table')))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite("data.table")
	library("data.table")
}
load('noodles.C.7.spaghetti.tumors.read.coverage.Rda')
tumor.read.stats.frame<-t(apply(noodles.C.7.spaghetti.tumors.read.coverage,1,quantile))
save(file='noodles.C.7.spaghetti.tumors.read.quantiles.Rda',list=c('tumor.read.stats.frame'))
