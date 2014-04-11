if (!suppressWarnings(require('rtracklayer')))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite("rtracklayer")
}

if (!suppressWarnings(require('DASiR')))
{
	source("http://bioconductor.org/biocLite.R")
	biocLite("DASiR")
}

vistaEnhancers.with.methylation.loaded<-FALSE
# we can the whole thing to vistaEnhancers.with.methylation.Rda
if(file.exists('vistaEnhancers.with.methylation.Rda'))
	if ('vistaEnhancers.with.methylation' %in% load('vistaEnhancers.with.methylation.Rda'))
		if (class(vistaEnhancers.with.methylation)=='data.frame')
			vistaEnhancers.with.methylation.loaded<-TRUE

if (!vistaEnhancers.with.methylation.loaded)
{
	source('../common/read_clinical.R')
	#Clinical prepared.

	#it is folder with bed files
	peakbedsfolder<-paste0(dataFolder,'/PeakCalls/bedfiles/')

	beds<-list.files(peakbedsfolder)

	#reading vistaEnhancers
	# we can the whole thing to vistaEnhancers.with.methylation.Rda
	source('../common/load_or_read_vista_enhancers.R')
	#print(vistaEnhancers)
	#vistaEnhancers read fro DAS or loaded
	vistaEnhancers.with.methylation<-as(vistaEnhancers,"data.frame")

	bed_available<-logical(0)
	bed_used<-rep(FALSE,length(beds))

	for (DNAid in DNAids)
	{
		DNAidKey<-strsplit(DNAid,',')[[1]][1]	#remove all after ,	
		match<-grep(DNAidKey,beds)
		if (!length(match)) 
		{
			DNAidKey<-paste0(strsplit(DNAid,'_')[[1]],collapse='') 
			#remove _ from key; sometimes, it help
			match<-grep(DNAidKey,beds)
		}
		if (!length(match)) 
		{
			bed_available<-c(bed_available,FALSE)
			next
		}
		message(DNAidKey)
		if (length(match)>1) stop(paste0("More than one match of DNAid ",DNAid," amonng the bed file names.\n"));
		bedfilename<-paste0(peakbedsfolder,beds[match[1]]);
		methylated.ranges<-as(import(bedfilename),"RangedData")
		overlaps<-findOverlaps(vistaEnhancers,methylated.ranges)
		methylcoverage=integer(0)
		for(chr in names(vistaEnhancers))
		#cycle by chromosome
		{
			list.of.ovelaps.in.this.chr<-as.list(overlaps[[chr]])
			width.of.meth.ranges.in.this.chr<-width(methylated.ranges[chr])
			methylcoverage.this.chr<-sapply(1:length(vistaEnhancers[chr][[1]]),function(band){
				sum(width.of.meth.ranges.in.this.chr[list.of.ovelaps.in.this.chr[[band]]])
			})#list of methylated coverage per cytoband
			methylcoverage<-c(methylcoverage,methylcoverage.this.chr)
			#add to main list
			#we need dual cycle because of the findOverlap return structure
		}
		message('done\n')
		vistaEnhancers.with.methylation[[DNAid]]=methylcoverage
		bed_used[match[1]]<-TRUE
		bed_available<-c(bed_available,TRUE)
	}
	message('Saving...')
	save(file='vistaEnhancers.with.methylation.Rda',list=c('vistaEnhancers.with.methylation','Clinical','clinFile','beds','bed_available','bed_used','tumors','normals','DNAids'))
}

tests.number<-dim(vistaEnhancers.with.methylation)[1]

wilcoxon.p.values<-numeric(tests.number)
normals.are.less.methylated<-logical(tests.number)

expected.w.statistic<-(sum(normals[bed_available])*sum(tumors[bed_available]))/2

for (rown in 1:tests.number)
{
	meth.values<-jitter(as.numeric(vistaEnhancers.with.methylation[rown,][DNAids[bed_available]]))
	w<-wilcox.test(meth.values[normals[bed_available]],meth.values[tumors[bed_available]])
	wilcoxon.p.values[rown]<-w$p.value
	normals.are.less.methylated[rown]<-(w[['statistic']]<expected.w.statistic)
	#anova.result<-c(anova.result,anova(lm(meth.values~tumors))[1,'Pr(>F)'])
}

wilcoxon.p.values.bonferroni<-p.adjust(wilcoxon.p.values,'bonferroni')
#wilcoxon.p.values.fdr<-p.adjust(wilcoxon.p.values,'fdr')

DM.enhancers<-which(wilcoxon.p.values<=0.05)
DM.enhancers.bonferroni<-which(wilcoxon.p.values.bonferroni<=0.05)
#DM.enhancers.fdr<-which(wilcoxon.p.values.fdr<=0.05)
#here, we form output statictics
columns<-c('space','start','end','id','score')
vistaEnhancers.stat<-cbind(vistaEnhancers.with.methylation[,columns],'p.value'=wilcoxon.p.values,'is.hyper'=normals.are.less.methylated)
DM.vistaEnhancers.stat<-vistaEnhancers.stat[DM.enhancers.bonferroni,]
message('Looking for closest genes')

source('../common/load_or_read_refseq_genes.R')

DM.vistaEnhancers.Ranges<-as(DM.vistaEnhancers.stat[,columns],'RangedData')

downstream<-character(0)
upstream<-character(0)

for (chr in names(DM.vistaEnhancers.Ranges))
{
	upstream.indices<-precede(DM.vistaEnhancers.Ranges[chr]$ranges,refseqGenes[chr]$ranges)
	downstream.indices<-follow(DM.vistaEnhancers.Ranges[chr]$ranges,refseqGenes[chr]$ranges)
	upstream<-c(upstream,as.character(refseqGenes[chr]$label)[upstream.indices])
	downstream<-c(downstream,as.character(refseqGenes[chr]$label)[downstream.indices])
}
message('done\n')
message('Mapping RefSeq to HGNC name')
source('../common/load_or_read_HGNC_ids.R')

HGNC_u_coord_list<-as.integer(sapply(upstream,
				function(RSid)
				{
					coord<-which(hgnc.ids$RefSeq.IDs==RSid)
					#id<-hgnc.ids$HGNC.ID[coord[1]]
					coord
				},USE.NAMES=FALSE
			))

HGNC_d_coord_list<-as.integer(sapply(downstream,
				function(RSid)
				{
					coord<-which(hgnc.ids$RefSeq.IDs==RSid)
					#id<-hgnc.ids$HGNC.ID[coord[1]]
					coord
				},USE.NAMES=FALSE
			))
message('done\n')
DM.vistaEnhancers.stat<-cbind(DM.vistaEnhancers.stat,'upsream'=upstream,'HGNC'=hgnc.ids$HGNC.ID[HGNC_u_coord_list],'name'=hgnc.ids$Approved.Symbol[HGNC_u_coord_list],'downstream'=downstream,'HGNC'=hgnc.ids$HGNC.ID[HGNC_d_coord_list],'name'=hgnc.ids$Approved.Symbol[HGNC_d_coord_list])

