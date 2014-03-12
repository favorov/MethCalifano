source('Y_chrom_methylation.R')

print(wilcox.test(Y.chrom.stats$Peak.length.sum[Y.chrom.stats$male],Y.chrom.stats$Peak.length.sum[!Y.chrom.stats$male]))
print(wilcox.test(Y.chrom.stats$Peak.score.sum[Y.chrom.stats$male],Y.chrom.stats$Peak.score.sum[!Y.chrom.stats$male]))
print(wilcox.test(Y.chrom.stats$Product.sum[Y.chrom.stats$male],Y.chrom.stats$Product.sum[!Y.chrom.stats$male]))
print(wilcox.test(one_m_j,zero_f_j))

