#the script work without parallel frame if you just run it
#worker:  worker my-worker-no workers-no (good for any parallel enviroment, including fork)
#combine: combiner workers-no (good for any parallel enviroment, including fork)
#everuthing happens in current directory!

noodle.lenght<-100
noodle.code<-'C'

args <- commandArgs(trailingOnly = TRUE)

i.am.worker<-FALSE
i.am.combiner<-FALSE
i.am.alone<-FALSE

if (length(args) > 0)
{
	if (! args[1] %in% c('worker','combiner'))
		stop('First argument is unknown!')
	if ('worker' == args[1] )
	{
		if (length(args)!=3)
			stop('worker is to have two more args')
		if(suppressWarnings(is.na(my.worker.no<-as.integer(args[2]))))
			stop('Second arg is not number')
		if(suppressWarnings(is.na(workers.no<-as.integer(args[3]))))
			stop('Third arg is not number')
		if( my.worker.no>500)
			stop('My number of worker is a strange (>500) number')
		if(workers.no<1)
			stop('Number of workers is a strange (<1) number')
		if(my.worker.no>workers.no)
			stop('Number of workers is less than my worker no')
		i.am.worker<-TRUE
	} else if ('combiner'==args[1])
	{
		if (length(args)!=2)
			stop('combiner is to have one more args')
		if(suppressWarnings(is.na(workers.no<-as.integer(args[2]))))
			stop('Second arg is not number')
		if( workers.no>500)
			stop('My number of worker is a strange (>500) number')
		if(workers.no<1)
			stop('Number of workers is a strange (<1) number')
		i.am.combiner<-TRUE
	} else
		stop('First arg is not a \'worker\' or \'combiner\'')
} else
	i.am.alone<-TRUE

if (length(which(c(i.am.alone,i.am.worker,i.am.combiner))) != 1)
	stop ('Something wrong with the self-identification of the script')

if(i.am.alone || i.am.combiner)
{	
	resultfilename<-paste0('noodles.',noodle.code,'.fisher.results.Rda')
	fisher.results.var.name<-paste0('fisher.noodles.',noodle.code,'.result')
} else #i.am.worker
{
	resultfilename<-paste0('noodles.',noodle.code,'.fisher.results.worker.',my.worker.no,'.Rda')
	fisher.results.var.name<-'fisher.noodles.result.mat'
}

noodles.fisher.results.loaded<-FALSE
# we can load the whole thing


if(file.exists(resultfilename))
	if (fisher.results.var.name %in% load(resultfilename))
			noodles.fisher.results.loaded<-TRUE
#if we loaded it, we do nothing.

if(!noodles.fisher.results.loaded)
{
	if(i.am.alone || i.am.worker) 
	{
		
		if (!suppressWarnings(require('differential.coverage')))
		{
			if (!suppressWarnings(require('devtools')))
			{
				source("http://bioconductor.org/biocLite.R")
				biocLite("devtools")
				library("devtools")
			}
			install_github('favorov/differential.coverage')
			#load_all('../../../../../differential.coverage/')
			library('differential.coverage')
		}

		#we need this to load - we are going to Fisherise
		noodles.loaded<-FALSE
		# we can load the whole thing from noodles.$noodle.code.Rda

		noodles.file<-paste0('noodles.',noodle.code,'.Rda')
		noodles.methylation.var.name<-paste0('noodles.',noodle.code,'.methylation') # it is a variable name!
		noodles.ranges.var.name<-paste0('noodles.',noodle.code)

		if(file.exists(noodles.file))
		{
			loaded<-load(noodles.file)
			if (noodles.methylation.var.name %in% loaded) 
				if (class(get(noodles.methylation.var.name))=='dgCMatrix')
					if (noodles.ranges.var.name %in% loaded)
						if(class(get(noodles.ranges.var.name))=='GRanges')
					noodles.loaded<-TRUE
		}
		if(!noodles.loaded)
		{
			if (i.am.worker) 
				stop (paste0('A worker (',my.worker.no,') cannot load noodles+methylation info from ',noodles.methylation.var.name,' and ',noodles.ranges.var.name,'.\n'))
			#if here, i.am.alone
			source(paste0('prepare.gw.noodles.',noodle.code,'.R'))
		}
		
		#assign(noodles.methylation.var.name,get(noodles.methylation.var.name)[1:60000,]) #test
		message('fishering')

		noodles.number<-dim(get(noodles.methylation.var.name))[1]
		#project-specific code; the bed.ids were prepared earlier
		#here, we suppose that contrast is prepare earlier, and loaded
		#and it is 0 for normal, 1 for tumor
		norm.no<-length(which(contrast==0))
		tumor.no<-length(which(contrast==1))

		fishtabs<-as.matrix(prepare.tabulated.fisher(tumor.no,norm.no))

		if (i.am.alone)
		{
			tests.number<-noodles.number
			my.worker.start=1
			my.worker.end=noodles.number
		}
		else
		{
			step.tests.number<-noodles.number %/% workers.no + ifelse(noodles.number %% workers.no > 0,1,0) #if remainder is zero, / is ok 
			my.worker.start<-1+step.tests.number*(my.worker.no-1)
			my.worker.end<-min(my.worker.start+step.tests.number-1,noodles.number)
			tests.number<-my.worker.end-my.worker.start+1 
			# real number of tests, critical for the last worker when my.worker.end==noodles.number 
			# and thus tests.number < step.tests.number

		}
	
		message('create result matrix')
		fisher.noodles.result.mat<-matrix(fishtabs[1,],ncol=6,nrow=tests.number,byrow=TRUE)
		
		colnames(fisher.noodles.result.mat)<-c('fisher.p.values','meth.in.normals.ratio','meth.in.tumors.ratio','OR','CI_95_L','CI_95_H') 
		
		revcontrast<-!contrast
		report.every<-tests.number %/% 100
		message('fill result matrix')

		for (rown in 1:tests.number) 	
		{
			therown<-rown+my.worker.start-1
			if (!(rown %% report.every)) message(rown)
			metraw<-get(noodles.methylation.var.name)[therown,]
			aslogic<-as.logical(metraw)
			MY<-sum(aslogic & contrast)
			MN<-sum(aslogic & revcontrast)
			if (0==MN && 0==MY) next
			fishres<-fishtabs[tab.fisher.row.no(tumor.no,norm.no,MY,MN),]
			fisher.noodles.result.mat[rown,]<-fishres
		}
		if (i.am.alone)
		{
			message('converting to dataframe')
			assign(fisher.results.var.name,as.data.frame(fisher.noodles.result.mat))
			message('done\n')
			message('Saving...\n')
			save(file=resultfilename,list=c(fisher.results.var.name,'tests.number','contrast'))
		} else #worker
		{
			message('Saving...\n')
			save(file=resultfilename,list=c(fisher.results.var.name,'noodles.number','my.worker.no','my.worker.start','workers.no','my.worker.end','contrast'))
		}
		message('done...\n')
		#in the worker case, 'fisher.noodles.result.mat' == fisher.results.var.name
	}
	else #combiner
	{
		message('Combiner started...\n')
		#testing the folder
	 	rdalist<-dir(pattern=paste0('noodles.',noodle.code,'.fisher.results.worker.[0-9]*.Rda'))
		if(length(rdalist)!=workers.no)
			stop(paste0('combiner: folder has other noodles.',noodle.code,'.fisher.results.worker.NN.Rda files (',length(rdalist),') than the workers.no (',workers.no,').'))
		worker.no<-1
		loadfilename<-paste0('noodles.',noodle.code,'.fisher.results.worker.',worker.no,'.Rda')
		load(loadfilename)
		prev.contrast<-contrast
		if (1 != my.worker.no)
			stop('Combiner error: first file has non-1 my.worker.no')
		if (1 != my.worker.start)
			stop('Combiner error: first file has non-1 my.worker.start')
		prev.worker.no<-my.worker.no
		prev.worker.end<-my.worker.end
		first.noodles.number<-noodles.number
		combinedresult<-fisher.noodles.result.mat
		for (worker.no in 2:workers.no)
		{
			loadfilename<-paste0('noodles.',noodle.code,'.fisher.results.worker.',worker.no,'.Rda')
			load(loadfilename)
			prev.contrast<-contrast
			if (worker.no != my.worker.no)
				stop('Combiner error: worker.no and number of file differ')
			if (prev.worker.end+1 != my.worker.start)
				stop('Combiner error:  my.worker.start is not prev.worker.end+1')
			if(first.noodles.number!=noodles.number)
				stop('Combiner error: noodles.number varies')
			prev.worker.no<-my.worker.no
			prev.worker.end<-my.worker.end
			combinedresult<-rbind(combinedresult,fisher.noodles.result.mat)
		}
		if (prev.worker.end != noodles.number)
			stop(paste0('Combiner error: the combined crowd is ', prev.worker.end, ' length while there are ', noodles.number,' noodles.'))
		message('converting to dataframe')
		assign(fisher.results.var.name,as.data.frame(combinedresult))
		message('done\n')
		message('Saving...\n')
		tests.number<-noodles.number
		save(file=resultfilename,list=c(fisher.results.var.name,'tests.number','contrast'))
	}
}
