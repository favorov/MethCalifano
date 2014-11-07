#prepare bed file names for Joe Califano's project
#it is folder with bed files
peakbedsfolder<-paste0(data.folder,'/PeakCalls/bedfiles/')

beds<-list.files(peakbedsfolder)

bed.used<-logical(length(beds))

DNAids<-Clinical$'Tumor/Tissue DNA HAND ID'[normals | tumors]
codes<-rownames(Clinical)[normals | tumors]
#we do not want to work with xeongraft & cell lines


contrast<-integer(0)

contrast[normals & ! tumors] <- 0
contrast[tumors & ! normals] <- 1
#contrast : 0 for normals, 1 for tumors, NA for unknown 

clinical.row.used<-!is.na(contrast)

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
		if (!length(match)) stop(paste0("No bed file name matched DNAid ",DNAid,".\n")); 
		if (length(match)>1) stop(paste0("More than one bed file name matches DNAid ",DNAid,".\n"));
		bedfilename<-paste0(peakbedsfolder,beds[match[1]]);
		bed..used<-get('bed.used',pos=1)
		bed..used[match[1]]=TRUE
		assign(x='bed.used',value=bed..used,pos=1)
		bedfilename
	}
)
message('Files assigned.\n')
bedfilenames<-unlist(bedfilenames)
