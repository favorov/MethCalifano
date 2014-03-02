source("Y_chrom_stats.R")
print(Y.chrom.stats)

print(wilcox.test(Y.chrom.stats$Peak.length.sum[Y.chrom.stats$male],Y.chrom.stats$Peak.length.sum[!Y.chrom.stats$male]))
print(wilcox.test(Y.chrom.stats$Peak.score.sum[Y.chrom.stats$male],Y.chrom.stats$Peak.score.sum[!Y.chrom.stats$male]))
print(wilcox.test(Y.chrom.stats$Product.sum[Y.chrom.stats$male],Y.chrom.stats$Product.sum[!Y.chrom.stats$male]))
#to comare: 
one_m<-rep(1,sum(Y.chrom.stats$male))
zero_f<-rep(0,sum(!Y.chrom.stats$male))
#remove ties by adding noise
one_m<-one_m+rnorm(length(one_m),0,0.00001)
zero_f<-zero_f+rnorm(length(zero_f),0,0.00001)
print(wilcox.test(one_m,zero_f))
