######### Final Script for ICD Regression Model #########

#In R studio

#set wd
setwd("/rdsgpfs/general")

#load fread (read in) file associated packages
library(data.table)
library(bit64)
library(PheWAS)

#### Step 1: Read in Genotype file, re-format Phenotype file and select out 'Sex'####

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

#read in phenotype with new type1 R data classes 
	#QC: Identify and convert 'fake' NA values into real NA values that can be detected by R using na.strings
d<-fread("/rdsgpfs/general/project/lms-ware-raw/live/UKBB/r/ukb29171.tab", header=TRUE, sep="\t", na.strings=c("Na","NA"), colClasses=type1) #3 mins
#select, bind and rename sample ID and 'Sex' columns and assign to 'sex'
sex<-as.data.frame(cbind(d$f.eid, d$f.31.0.0),stringsAsFactors=FALSE)
colnames(sex)<-c("id","sex")
sex$id<-as.integer(sex$id)
levels(sex[,2])<-c("M","F")

#read in new genotype file for 'additive' genetic model instead of dominant
c <- read.table("/rdsgpfs/general/project/lms-ware-analysis/live/xiao/M.sc/genotype/genotype_add.txt", header = TRUE, sep = " ", na.strings=c("Na","NA"))
#assign Genotype for PheWAS package 
genotypes=as.data.frame(e[,c(1,3)],stringsAsFactors=FALSE) #column 2 for BAG3, column 3 for FHOD3
colnames(genotypes)[1]<-"id"
genotypes$id<-as.integer(genotypes$id)

#### Step 2: Load UKBB ICD file and select relavent columns ####

#Load UKBB additional ICD Phenotype file raw data and filter out relavent samples
f<-fread("/rdsgpfs/general/project/lms-ware-raw/live/UKBB/csv/additional_1/ukb37190.csv", header=TRUE, sep=",", na.strings=c("Na","NA"), colClasses="character")
#d1 <- d[order(eid),] 

	#QC: Automatic filter to remove background non-caucasian population by ethnicity
	pop<-cbind(d$f.eid, d$f.22006.0.0) #get ethnicity information from UKBB phenotype file, 1 is Caucasian. NA others
	colnames(pop)<-c("eid","pop") #rename column variables
	f1 <- merge(f, pop, by = "eid", all = TRUE, sort = FALSE) #merge ethinicity information to UKBB ICD file by sample ID
	f1<-f1[!is.na(f1$pop), ] #select out unique values i.e exclude rows with NA values

#### Step 3: select relevant columns and Re-format ICD codes from the selected columns from the UKBB file into correct ICD-10-CM codes ####

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

#### Step 3: Regression analysis between BAG3 genotype data and the re-formatted ICD phenotype data using PheWAS package ####


#Re-shape and define phenotype functions to a 'wide' logical format
phenotypes=createPhenotypes(s1, id.sex=sex, min.code.count=1,vocabulary.map=PheWAS::phecode_map_icd10,,translate=TRUE)

#phenotypes=reshape(id.vocab.code.count[,-2],idvar="id", timevar="code",direction="wide")
#phenotypes[phenotypes=="20"]<-"TRUE"
#phenotypes[is.na(phenotypes)==TRUE]<-"FALSE"

	#QC: add non-record samples
	id_p<-phenotypes$id #phenotype ID
	id_a<-f1$eid #all sample IDs
	id_d<-setdiff(id_a,id_p) #see difference in ID, remove duplicates
	add<-as.data.frame(matrix(FALSE,length(id_d),ncol(phenotypes)))
	add[,1]<-id_d
	colnames(add)<-colnames(phenotypes)
	phenotypes<-rbind(phenotypes,add)
	phenotypes[,1]<-as.integer(phenotypes[,1])

#output intermediate phenotype file
write.table(phenotypes,"/rdsgpfs/general/project/lms-ware-analysis/live/xiao/M.sc/ICD_wide.txt",row.names=FALSE,quote=FALSE,col.names=T,sep="\t") #output table
phenotypes<-fread("/rdsgpfs/general/project/lms-ware-analysis/live/xiao/M.sc/ICD_wide.txt",header=TRUE,sep="\t") #read in data summary file


#define complete combined dataset
data=inner_join(inner_join(phenotypes,genotypes),sex)

#run PheWAS package and get results 
results=phewas_ext(data,phenotypes=names(phenotypes),genotypes=c("FHOD3"),covariates=c("sex"), cores=1,min.records=20)

#output data table of results
write.table(results,"/rdsgpfs/general/project/lms-ware-analysis/live/xiao/M.sc/ICD_result_FHOD3.txt",row.names=FALSE,quote=FALSE,col.names=TRUE,sep="\t") 
results<-read.table("/rdsgpfs/general/project/lms-ware-analysis/live/xiao/M.sc/ICD_result_FHOD3.txt",header=TRUE,sep="\t") #read in data summary file


#output plot of results
png("/rdsgpfs/general/project/lms-ware-analysis/live/Prarthna/ICD_result_FHOD3plot.png", res=150, width=1000, height=1000)
phewasManhattan(results, title="FHOD3", annotate.size=3, max.x=1500)
dev.off()