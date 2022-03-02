##################################################################################################################"
##################                                                                         #######################"
##################           Preparing genomic data for the DRYAD repository               #######################"
##################                                                                         #######################"
##################################################################################################################"


# Here I format the genomic data as in the reports:
#  - GenomicRelationshipMatrices.Rmd
#  - CalculatingGPEAandRPEA.Rmd


# Loading the file with the genotype of each clone for each SNP 
# => 5165 SNPs in rows and 523 genotypes (i.e. clones) in columns
geno <- read.csv("data/5165snps523genotypesNA.txt", header=FALSE, row.names=1) 


# Removing the first two columns with allele info (A,T, G or C)
geno <- geno[,3:dim(geno)[[2]]]


# In this file, SNPs have their names, but not the genotypes, that's why we also load a file with the genotype names (clone names)
geno_names <- read.delim2("data/ClonapinBlups523IndPiMASSJuly2019.txt", row.names=1)


# We give the genotype name for each column of geno
colnames(geno) <- rownames(geno_names)


# We save the genomic data and this dataset will be the one in the DRYAD repository
write.csv(geno,
          file= paste0("data_DRYAD/GenomicData_",nrow(geno),"SNPs_",ncol(geno),"clones.csv"),
          row.names = T)
