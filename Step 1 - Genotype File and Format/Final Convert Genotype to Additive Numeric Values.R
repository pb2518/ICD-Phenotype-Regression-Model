######## Function: Convert genotype to numeric format based on genetic model ########

#### Parameter of genetic model ####

#Can be additive, dominant, recessive
m<-"additive"
#m<-"dominant"
#m<-"recessive"

#### Step 1: Load inputs ####

#snp_list file with risk allele information
snp<-as.matrix(read.table("/rdsgpfs/general/project/lms-ware-analysis/live/xiao/M.sc/genotype/geno_info.txt",header=TRUE))

#sample genotype file
geno<-as.matrix(read.table("/data/Share/xiao/HCM_GWAS/analysis/3D_phenotype/Msc/HCM_snp_merged.ped",header=FALSE))
geno1<-geno[,7:ncol(geno)] #select only genotype information

#### Step 2: Get the number of alleles for each snp, for each individual, e.g: 0,1 ####

a<-matrix(,nrow(geno),nrow(snp)) #define empty matrix for output
for (i in 1:nrow(snp)) #loop for snp
{
allele_ref<-snp[i,"risk_allele"] #get risk allele
	for (j in 1:nrow(geno)) #loop for sample
	{
	allele<-c(geno1[j,2*i-1],geno1[j,2*i]) #get allele for each snp
		if (("0" %in% allele)==FALSE){ #if no missing genotype
			if (m=="additive"){a[j,i]<-length(which(allele==allele_ref))}
			if (m=="dominant"){a[j,i]<-a[j,i]<-as.numeric(allele_ref %in% allele)}
			if (m=="recessive"){a[j,i]<-as.numeric(length(which(allele==allele_ref))==2)}
		} else {
		a[j,i]<-NA #if presence of missing genotype
		}
	
	}
}

#out<-cbind(geno[,1:6],a)
out<-cbind(geno[,1],a)
colnames(out)<-c("f.eid","BAG3","FHOD3")

#### Step 3: Output table ####
write.table(out,"/rdsgpfs/general/project/lms-ware-analysis/live/xiao/M.sc/genotype/genotype_add.txt",col.names=TRUE, row.names=FALSE,quote=FALSE)
