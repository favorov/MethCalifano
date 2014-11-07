#reads clinical file for Joe Califano's projects

if ((! 'data.version' %in% ls()) || is.na(data.version)) data.version=5


if (data.version!=4 && data.version!=5 && data.version != 'final') stop(paste0('Unknown data version',data.version))


if ((! 'data.folder' %in% ls()) || is.na(data.folder)) data.folder<-'../../Data'

#the file to be read
if(data.version==4)	
{
	clinFile <- paste0(data.folder,'/Clinical/2014-02-06_RO1_Batch1-4_ClinicalData_DeID.csv')
	Clinical <- read.csv(clinFile,stringsAsFactors=F)

	#fix up the column names to account for the fact that there are two header rows
	#the names are absent for first three rows, they were in a row below 
	colnames(Clinical)[1:3] <- Clinical[1,1:3]
	Clinical <- Clinical[2:nrow(Clinical),]

	#remove all service rows (where the TissueDNAID is "")
	Clinical<-Clinical[!Clinical[,2]=="",]

	#rownames are Codes
	rownames(Clinical)<-Clinical$Codes

	# determine which are tumor and which are normal samples
	tumors <- as.integer(Clinical$Code) < 99
	normals <- as.integer(Clinical$Code) > 99

	DNAids<-Clinical$'Tissue DNA HAND ID'
	#clinical prepared
}

if(data.version==5 || data.version=='final')	
{
	clinFile <- paste0(data.folder,'/Clinical/2014-10-07_DiscoveryCohort_ClinicalData-Final_DeID.csv')
	Clinical <- read.csv(clinFile,stringsAsFactors=F)

	#fix up the column names to account for the fact that there are two header rows
	#the names are absent for first three rows, they were in a row below 
	colnames(Clinical)[1:3] <- Clinical[1,1:3]
	Clinical <- Clinical[2:nrow(Clinical),]

	#remove all service rows (where the TissueDNAID is "")
	Clinical<-Clinical[!Clinical[,2]=="",]

	#rownames are Codes
	rownames(Clinical)<-Clinical$code

	# determine which are tumor and which are normal samples
	tumors <- as.integer(Clinical$code) < 99
	normals <- (as.integer(Clinical$code) > 99) & ((as.integer(Clinical$code) < 200) )

	#clinical prepared
}

