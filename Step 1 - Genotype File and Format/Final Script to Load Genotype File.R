######## Function: Re-format and Load genotype table and output required columns for analysis #########

#In R studio

#### Step 1: Load Genotype file for 'additive' model ####

#swd
setwd("/rdsgpfs/general")

#load genotype file (numeric 'additive' genotype version, see script used to convert) 
	#QC: Identify and convert 'fake' NA values into real NA values that can be detected by R using na.strings
a <- read.table("/rdsgpfs/general/project/lms-ware-analysis/live/xiao/M.sc/genotype/genotype_add.txt", header = TRUE, sep = " ", na.strings=c("Na","NA"))

	#QC: check structure and format
	head(a, n = 10)
	str(a) 
	#QC: check data class
	class(a$f.eid)
	class(a$BAG3)
	class(a$FHOD3)


#For BAG3 (chr10,rs2234962 ), 0 = TT, 1 = CT/TC, and 2 = CC (Additive for C)
#For FODH3 (chr18, rs2303510), 0 = GG, 1 = AG/GA, and 2 = AA (Additive for A)


#### Step 2: Select out and assign genotype for analysis ####

#assign Genotype for PheWAS package 
genotypes=as.data.frame(e[,c(1,2)],stringsAsFactors=FALSE) #column 2 for BAG3, or column 3 for FHOD3
colnames(genotypes)[1]<-"id" #re-name sample ID column
genotypes$id<-as.integer(genotypes$id) #Re-format datatype