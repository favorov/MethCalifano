library(data.table)
load('noodles.C.7.spaghetti.tumors.read.coverage.Rda')
tumor.read.stats.frame<-t(apply(noodles.C.7.spaghetti.tumors.read.coverage,1,quantile))
save(file='tumor.read.noodle.stats.Rda',list=c('tumor.read.stats.frame'))
