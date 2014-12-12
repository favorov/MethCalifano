#we suppose all the preparatio is already done
library('Rsamtools')
library('GenomicAlignments')
load('noodles.C.Rda')
BAM.folder<-'f:/Califano/MBDseq_PeakCalls/bam-bai_files/'
b <- readGAlignments('f:/Califano/MBDseq_PeakCalls/bam-bai_files/002_38623_Enrich_hg19.combined.bam');
over<-countOverlaps(noodles.C,b)


