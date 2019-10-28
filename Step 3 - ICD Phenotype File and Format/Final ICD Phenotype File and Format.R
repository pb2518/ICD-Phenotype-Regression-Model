######## Function: Re-format ICD Phenotype codes and load into a table #########

#In R studio

#set wd
setwd("/rdsgpfs/general")

#load fread (read in) file associated packages
library(data.table)
library(bit64)

#### Step 1: Load raw UKBB ICD Phenotype file and filter Out relavent columns ####

#Load UKBB additional ICD Phenotype file raw data and filter out relavent samples
d<-fread("/rdsgpfs/general/project/lms-ware-raw/live/UKBB/csv/additional_1/ukb37190.csv", header=TRUE, sep=",", na.strings=c("Na","NA"), colClasses "character")
#d1 <- d[order(eid),] 

	#QC: Automatic filter to remove background non-caucasian population by ethnicity
	pop<-cbind(d$f.eid, d$f.22006.0.0) #get ethnicity information from UKBB phenotype file, 1 is Caucasian. NA others
	colnames(pop)<-c("eid","pop") #rename column variables
	f1 <- merge(f, pop, by = "eid", all = TRUE, sort = FALSE) #merge ethinicity information to UKBB ICD file by sample ID
	f1<-f1[!is.na(f1$pop), ] #select out unique values i.e exclude rows with NA values

#### Step 2: select relevant columns and Re-format ICD codes from the selected columns from the UKBB file into correct ICD-10-CM codes ####

#select first 214 columns of diagnosis data
ICD<-as.matrix(f1[350001:nrow(f1),1:214]) #run script in sections to reduce run time as follows (takes 1 hr +)
#p1,1:50000
#p2,50001:100000
#p3,100001:150000
#p4,150001:200000
#p5,200001:250000
#p6,250001:300000
#p7,300001:350000
#p8,350001:nrow(f1)


#create a function that identifies ICD code from UKBB ICD file
getICD <- function(x) {
  a<-unique(x)
  a<-a[!a == ""]
  a1<-a[1]
  a<-a[!a==a[1]]
  b<-cbind(a1,"ICD10CM",a,"20") #value 20 given to index 
  if (length(a)!=0) {return(b)}
}

s1<-c("id","vocabulary_id","code","index") #convert disease information into long format
for (i in 1:nrow(ICD)) #select out these ICD codes using the getICD function created using a loop for each sample
{
s<-getICD(ICD[i,])
s1<-rbind(s1,s)
print(i)
}

s1<-s1[-1,] #remove 1st header row

#output each table and use cat() to combine into one table
write.table(s1,"/rdsgpfs/general/project/lms-ware-analysis/live/xiao/M.sc/ICD_p8.txt",row.names=FALSE,quote=FALSE,col.names=F,sep="\t") 
s1<-fread("/rdsgpfs/general/project/lms-ware-analysis/live/xiao/M.sc/ICD_all.txt",header=F,sep="\t") #read in data summary file

colnames(s1)<-c("id","vocabulary_id","code","index") #assign column names  

#create translation table： UKBB code to ICD-10-CM code, using coding19 file downloaded from UKBB website
code <- as.data.frame(fread("/rdsgpfs/general/project/lms-ware-analysis/live/xiao/M.sc/coding19.tsv", header = TRUE, sep = "\t",stringsAsFactors=FALSE))
cc1<-strsplit(code[,2]," ")
cc2<-unlist(lapply(cc1, `[[`, 1))
cc3<-cbind(code[,1],cc2)

#translate all the UKBB ICD codes that do not match their corresponding ICD-10-CM code into the correct ICD-10-CM code
s1[,3]<-cc3[match(s1[,3], cc3[,1], incomparables = NULL),2]

#correct column datatypes
s1=as.data.frame(s1,,stringsAsFactors=FALSE)
s1$id<-as.integer(id.vocab.code.count$id)
s1$index<-as.numeric(id.vocab.code.count$index)
#s1[,2]<-"ICD10"
