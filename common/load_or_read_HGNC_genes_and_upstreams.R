#combines uscs gened and HGNC ids.
#writes tables for direc and reverse genes  and for upstreams 
#saves the information to ../common/CpGIs.Rda
#if the file exists, read it after some checks

source('../common/load_or_read_chrom_ranges.R')
source('../common/load_or_read_known_genes.R')
source('../common/load_or_read_HGNC_ids.R')
