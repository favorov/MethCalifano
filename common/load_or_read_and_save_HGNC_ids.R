#fetches HGNC ids  from HGNC server
#saves the information to ../common/CpGIs.Rda
#if the file exists, read it after some checks

HGNC.ids.loaded<-FALSE

if(file.exists('../common/HGNCids.Rda'))
	if ('hgnc.ids' %in% load('../common/HGNCids.Rda'))
		if (class(hgnc.ids)=='data.frame')
			HGNC.ids.loaded<-TRUE


if(!HGNC.ids.loaded)
{
	HGNC.url<- paste0('http://www.genenames.org/cgi-bin/download?',
          'col=gd_hgnc_id&',
          'col=gd_status&',
          'col=gd_locus_type&',
          'col=gd_pub_chrom_map&',
          'col=gd_pub_refseq_ids&',
          'col=md_ucsc_id&',
          'status=Approved&',
          'status=Entry%20Withdrawn&',
          'status_opt=2&',
          'where=&',
          'order_by=gd_hgnc_id&',
          'format=text&',
          'limit=&',
          'submit=submit')
	hgnc.ids<-read.table(HGNC.url,header=TRUE,sep="\t",stringsAsFactors=FALSE,fill=TRUE,quote = "\"")
	save(file='HGNCids.Rda',list=c('hgnc.ids'))
}

