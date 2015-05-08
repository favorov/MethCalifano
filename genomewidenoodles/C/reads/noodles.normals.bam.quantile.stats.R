#library(data.table)
load('noodles.C.7.spaghetti.normals.read.coverage.Rda')
norm.read.stats.frame<-t(apply(noodles.C.7.spaghetti.normals.read.coverage,1,quantile))
save(file='noodles.C.7.spaghetti.nolmals.read.quantiles.Rda',list=c('norm.read.stats.frame'))
