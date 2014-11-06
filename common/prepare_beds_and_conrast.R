#prepare bed file names for Joe Califano's project
#it is folder with bed files
peakbedsfolder<-paste0(data.folder,'/PeakCalls/bedfiles/')

beds<-list.files(peakbedsfolder)

DNAids<-DNAids[normals | tumors]
#we do not want to work with xeongraft & cell lines

contrast<-integer(0)

contrast[normals & ! tumors] <- 0
contrast[tumors & ! normals] <- 1
#contrast : 0 for normals, 1 for tumors, NA for unknown 

contrast<-contrast[!is.na(contrast)]
#we do not want to work with xeongraft & cell lines

bedfilenames<-lapply(DNAids,function(DNAid){ 
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
		#message(match[1])
		#message(bedfilename)
		bedfilename
	}
)
message('Files assigned.\n')
bedfilenames<-unlist(bedfilenames)
