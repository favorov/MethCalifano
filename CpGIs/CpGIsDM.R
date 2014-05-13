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

CpGIs.with.methylation.loaded<-FALSE
# we can the whole thing to CpGIs.with.methylation.Rda
if(file.exists('CpGIs.with.methylation.Rda'))
	if ('CpGIs.with.methylation' %in% load('CpGIs.with.methylation.Rda'))
		if (class(CpGIs.with.methylation)=='data.frame')
			CpGIs.with.methylation.loaded<-TRUE

if (!CpGIs.with.methylation.loaded)
{
	source('../common/read_clinical.R')
	#Clinical prepared.

	#it is folder with bed files
	peakbedsfolder<-paste0(dataFolder,'/PeakCalls/bedfiles/')

	beds<-list.files(peakbedsfolder)

	#reading islands 
	# we can the whole thing to CpGIs.with.methylation.Rda
	source('../common/load_or_read_CpGIs.R')
	#islands are read fro DAS or loaded

	CpGIs.with.methylation<-as(CpGIs,"data.frame")

	bed_available<-logical(0)
	bed_used<-rep(FALSE,length(beds))

	for (DNAid in DNAids)
	{ 
		DNAidKey<-strsplit(DNAid,',')[[1]][1]	#remove all after ,	
		message(DNAidKey)
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
		if (length(match)>1) stop(paste0("More than one match of DNAid ",DNAid," amonng the bed file names.\n"));
		bedfilename<-paste0(peakbedsfolder,beds[match[1]]);
		methylated.ranges<-as(import(bedfilename),"RangedData")
		overlaps<-findOverlaps(CpGIs,methylated.ranges)
		methylcoverage=integer(0)
		for(chr in names(CpGIs))
		#cycle by chromosome
		{
				list.of.ovelaps.in.this.chr<-as.list(overlaps[[chr]])
				width.of.meth.ranges.in.this.chr<-width(methylated.ranges[chr])
				methylcoverage.this.chr<-sapply(1:length(CpGIs[chr][[1]]),function(island_no){
					#attention: length(CpGIs[chr][[1]] supposes that there is column in datarange other than space and ranges (e.g., Id)
					#otherwise, use  length(start(CpGIs[chr]))
					sum(width.of.meth.ranges.in.this.chr[list.of.ovelaps.in.this.chr[[island_no]]])
				})#list of methylated coverage per cytoband
			#list of methylated coverage per island
			#the main idea is that as.list(hitsobject) convert is to list of vectors, 
			#list is addressd by query oblects and give ref object lists (possibly, empty)
			methylcoverage<-c(methylcoverage,methylcoverage.this.chr)
			#add to main list
			#we need dual cycle because of the findOverlap return structure
		}
		CpGIs.with.methylation[[DNAid]]=methylcoverage
		bed_used[match[1]]<-TRUE
		bed_available<-c(bed_available,TRUE)
		message('done\n')
	}
	message('Saving...\n')
	save(file='CpGIs.with.methylation.Rda',list=c('CpGIs.with.methylation','Clinical','clinFile','beds','bed_available','bed_used','tumors','normals','DNAids'))
}

CpGIs.wilcoxon.data.loaded<-FALSE
# we can the whole thing to CpGIs.with.methylation.Rda
if(file.exists('CpGIs.wilcoxon.data.Rda'))
	if ('wilcoxon.p.values' %in% load('CpGIs.wilcoxon.data.Rda'))
			CpGIs.wilcoxon.data.loaded<-TRUE
if(!CpGIs.wilcoxon.data.loaded)
{
	message('Wilcoxon\n')
	wilcoxon.p.values<-numeric(0)

	normals.are.less.methyl.covered<-logical(0)

	expected.w.statistic<-(sum(normals[bed_available])*sum(tumors[bed_available]))/2

	tests.number<-dim(CpGIs.with.methylation)[1]

	for (rown in 1:tests.number)
	{
		meth.values<-jitter(as.numeric(CpGIs.with.methylation[rown,][DNAids[bed_available]]))
		w<-wilcox.test(meth.values[normals[bed_available]],meth.values[tumors[bed_available]])
		wilcoxon.p.values<-c(wilcoxon.p.values,w$p.value)
		normals.are.less.methyl.covered<-c(normals.are.less.methyl.covered,(w[['statistic']]<expected.w.statistic))
		#anova.result<-c(anova.result,anova(lm(meth.values~tumors))[1,'Pr(>F)'])
	}

	message('done\n')
	message('Saving...\n')
	save(file='CpGIs.wilcoxon.data.Rda',list=c('wilcoxon.p.values','normals.are.less.methyl.covered','tests.number'))
}

CpGIs.fisher.data.loaded<-FALSE
# we can the whole thing to CpGIs.with.methylation.Rda
if(file.exists('CpGIs.fisher.data.Rda'))
	if ('fisher.p.values' %in% load('CpGIs.fisher.data.Rda'))
			CpGIs.fisher.data.loaded<-TRUE
if(!CpGIs.fisher.data.loaded)
{
	message('Fishering\n')
	tests.number<-dim(CpGIs.with.methylation)[1]
	fisher.p.values<-numeric(tests.number)
	meth.in.normals.ratio<-numeric(tests.number)
	meth.in.tumors.ratio<-numeric(tests.number)
	OR<-numeric(tests.number)
	CI_95_L<-numeric(tests.number)
	CI_95_H<-numeric(tests.number)


	for (rown in 1:tests.number)
	{
		cotable<-table(as.logical(as.numeric(CpGIs.with.methylation[rown,][DNAids[bed_available]])),tumors[bed_available])
		if(nrow(cotable)==1)#nonmeth
		{
			fisher.p.values[rown]<-1.
			meth.in.tumors.ratio[rown]<-0
			meth.in.normals.ratio[rown]<-0
			OR[rown]<-NA
			CI_95_L[rown]<-NA
			CI_95_H[rown]<-NA
			next
		}
		fisherres<-fisher.test(cotable)
		fisher.p.values[rown]<-fisherres$p.value
		meth.in.tumors.ratio[rown]<-cotable[2,2]/cotable[1,2]
		meth.in.normals.ratio[rown]<-cotable[2,1]/cotable[1,1]
		OR[rown]<-fisherres$estimate
		CI_95_L[rown]<-fisherres$conf.int[1]
		CI_95_H[rown]<-fisherres$conf.int[2]
	}

	message('done\n')
	message('Saving...\n')
	save(file='CpGIs.fisher.data.Rda',list=c('fisher.p.values','tests.number','meth.in.tumors.ratio','meth.in.normals.ratio','OR','CI_95_L','CI_95_H'))
}

DM.W.CpGIslands<-which(wilcoxon.p.values<=0.05)
DM.W.CpGIslands.Bonferroni<-which(wilcoxon.p.values*tests.number<=0.05)

DM.F.CpGIslands<-which(fisher.p.values<=0.05)
DM.F.CpGIslands.Bonferroni<-which(fisher.p.values*tests.number<=0.05)
#here, we form output statictics

#bonferroni
DM.CpGIslands.Bonferroni<-union(DM.W.CpGIslands.Bonferroni,DM.F.CpGIslands.Bonferroni)

columns<-c('id','space','start','end')
DM.CpGIs.stat<-cbind(
	CpGIs.with.methylation[DM.CpGIslands.Bonferroni,columns],
	'wilcoxon.p.value'=wilcoxon.p.values[DM.CpGIslands.Bonferroni],
	'hyper?'=normals.are.less.methyl.covered[DM.CpGIslands.Bonferroni],
	'fisher.p.value'=fisher.p.values[DM.CpGIslands.Bonferroni],
	'tmr.ratio'=meth.in.tumors.ratio[DM.CpGIslands.Bonferroni],
	'nor.ratio'=meth.in.normals.ratio[DM.CpGIslands.Bonferroni],
	'OR'=OR[DM.CpGIslands.Bonferroni],
	'CI_95_L'=CI_95_L[DM.CpGIslands.Bonferroni],
	'CI_95_H'=CI_95_H[DM.CpGIslands.Bonferroni]
)

rownames(DM.CpGIs.stat)<-NULL

#we want to put each diffmet CpGi to a cytoband
source('../common/load_or_read_karyotype.R')
DM.CpGIs.Ranges<-as(DM.CpGIs.stat[,columns],'RangedData')


message('Mapping to karyotype...\n')

CpGIs.to.karyotype<-findOverlaps(DM.CpGIs.Ranges,karyotype,type="within")
cytobands.of.DM.cpgis=character(0)
for (chr in names(CpGIs.to.karyotype))
{
  #chromosome cycle
	len<-length(DM.CpGIs.Ranges[chr][[1]])
	if (len==0) next
	#we did not do it in 'big' cycles over
	#overlaps, because their queries were like cytobands, etc - no empty chromosomes
	cytoband.numbers.this.chr<-sapply(1:len,function(island_no){
		my_cytoz<-as.list(CpGIs.to.karyotype[[chr]])[[island_no]]
		if (!length(my_cytoz)) {NA} else {my_cytoz[1]}
	})
	cytobands.of.DM.cpgis<-c(cytobands.of.DM.cpgis,as.character(karyotype[chr][[1]][cytoband.numbers.this.chr]))
}

popdir<-getwd()
setwd('../CytoBands/') # not to recalculate all the CytoBand caches
source('../CytoBands/CytoBandDM.R')
setwd(popdir)

DM.CpGIs.stat<-cbind(DM.CpGIs.stat,'cytoband'=cytobands.of.DM.cpgis,'DM.band?'=cytobands.of.DM.cpgis%in%as.character(karyotype.with.methylation$id)[DM.cytobands])
message('done\n')


message('Looking for closest gene')

source('../common/load_or_read_refseq_genes_with_HGNC_id.R')

nearestTSS<-character(0)
strand<-character(0)
position<-integer(0)
distance<-integer(0)

for (chr in names(DM.CpGIs.Ranges))
{
	nearest.HGNC.TSS.indices<-nearest(DM.CpGIs.Ranges[chr]$ranges,refseqTSSWithHGNCids[chr]$ranges)
	nearestTSS<-c(nearestTSS,as.character(refseqTSSWithHGNCids[chr]$symbol)[nearest.HGNC.TSS.indices])
	this_chr_positions<-start(refseqTSSWithHGNCids[chr]$ranges[nearest.HGNC.TSS.indices])
	this_chr_strand<-as.character(refseqTSSWithHGNCids[chr]$orientation)[nearest.HGNC.TSS.indices]
	strand<-c(strand,this_chr_strand)
	position<-c(position,this_chr_positions)
	this_chr_distances_u<-this_chr_positions-end(DM.CpGIs.Ranges[chr]$ranges)
	this_chr_distances_d<-this_chr_positions-start(DM.CpGIs.Ranges[chr]$ranges)
	this_chr_distances<-ifelse(this_chr_distances_u>0,this_chr_distances_u,ifelse(this_chr_distances_d<0,this_chr_distances_d,0))
	this_chr_distances<-as.integer(ifelse(this_chr_strand=='-',this_chr_distances,-this_chr_distances))
	#if thread is + and positive distance means CpGi down from the gene, so we change it to opposite
	distance<-c(distance,this_chr_distances)
}
message('done\n')

DM.CpGIs.stat<-cbind(DM.CpGIs.stat,'TSS near'=nearestTSS,'pos'=position,'strand'=strand,'distance'=distance)

DM.CpGIs.stat$id<-substr(DM.CpGIs.stat$id,6,1000) # 1000 'any'; we strip first 'CpGi: ' from the id

write.table(DM.CpGIs.stat,file='DM.CpGIs.stat.tsv',sep='\t',row.names=FALSE)

if(!require('xtable'))
{
  source("http://bioconductor.org/biocLite.R")
  biocLite("xtable")
  library("xtable")  
}

if(file.exists("DM.CpGIs.Ranges.html")) {file.remove("DM.CpGIs.Ranges.html")}

print(xtable(DM.CpGIs.stat,digits=c(0,0,0,0,0,8,0,8,2,2,2,2,2,0,0,0,0,0,0), display=c('d','s','s','d','d','g','s','g','f','f','f','f','f','s','s','s','d','s','d')), type="html", file="DM.CpGIs.stat.html")
