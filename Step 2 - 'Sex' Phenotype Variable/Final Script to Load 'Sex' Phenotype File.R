######## Function: Re-format and select 'Sex' from phenotype table and load output file #########

#In R studio

#set wd
setwd("/rdsgpfs/general")

#load fread (read in) file associated packages
library(data.table)
library(bit64)

#### Step 1: Reformat data type for each data class ####

#load and select the correct column list of current data types for each column of the UKBB phenotype data
type<-read.table("/rdsgpfs/general/user/pb2518/home/data_type.txt",header=TRUE,sep="\t") #read in data type summary file
type1<-as.character(type[,4]) #select out data type column information

#convert data types into R data classes
type1[type1 == "Sequence"]<-"character"
type1[type1 == "Categorical (single)"]<-"factor"
type1[type1 == "Categorical (multiple)"]<-"factor"
type1[type1 == "Text"]<-"character"
type1[type1 == "Integer"]<-"integer" 
type1[type1 == "Date"]<-"character"
type1[type1 == "Continuous"]<-"numeric"
type1[type1 == "Curve"]<-"character"

#### Step 2: Load phenotype file ####

#read in phenotype file with new type1 R data classes
	#QC: Identify and convert 'fake' NA values into real NA values that can be detected by R using na.strings
d <- fread("/rdsgpfs/general/project/lms-ware-raw/live/UKBB/r/ukb29171.tab", header=TRUE, sep="\t", na.strings=c("Na","NA"), colClasses=type1) #3 mins


	#QC: Check data classes coerced correctly by comparing specific columns in 'd' to original datatypes in 'type'
	head(type,n=3) #row 2 with UDI 31-0.0 should be categorical (single)
	class(d$f.31.0.0) #should now be factor
	#QC: Check structure and layout of first patient ID column and second column "sex" then a few more columns to see if correct
	head(d[ ,c(1,2)], n=10)
	head(d[ ,c(1:4)], n=10)
	
	#QC: Check data type and number of missing values
	d1<-as.data.frame(d)#convert to data frame
	f1t<-lapply(d1, class)#get data classes of d1
	
	check_t<-rep(0, ncol(f1))
	check_na<-rep(0, ncol(f1))
	for (m in 1:length(type1))
	{
	print(paste(m,":",f1t[m],"/",type1[m]))#print data type
	print(paste(type[m,3],"/",length(which(!is.na(f[,..m]))),sep=" ")) #print number of missing values
		if (f1t[m]==type1[m])#check data type match
		{
		check_t[m]<-1
		}    
		
		if (type[m,3]==length(which(!is.na(f[,..m]))))#check na vaule match
		{
		check_na[m]<-1
		}    
	}
	print(paste("there are",sum(check_t),"/",ncol(f1),"columns with correct data type"))
	print(paste("there are",sum(check_na),"/",ncol(f1),"columns with correct NA values")) 

#### Step 3: Select out, re-format and assign sample ID and Sex columns ####

sex<-as.data.frame(cbind(d$f.eid, d$f.31.0.0),stringsAsFactors=FALSE) #bind columns and assign 'sex' to a dataframe table
colnames(sex)<-c("id","sex") #re-name selected columns
sex$id<-as.integer(sex$id) #re-format data type of sample ID
levels(sex[,2])<-c("M","F") #assign 'sex' level identifiers

