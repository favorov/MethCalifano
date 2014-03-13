dataFolder<-'../../Data'

#the file to be read
#clinFile <- paste0(dataFolder,'/Clinical/2013-11-25_RO1_Batch1_ClinicalData_DeID.xls.csv')
clinFile <- paste0(dataFolder,'/Clinical/2014-02-06_RO1_Batch1-4_ClinicalData_DeID.csv')
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
