##################################################################################################################"
##################                                                                         #######################"
##################           Extracting rPEAs and gPEAs names                              #######################"
##################                                                                         #######################"
##################################################################################################################"


library(readr)   # CRAN v1.3.1
library(tidyverse)

# Load the genomic data:
geno <- read.csv("~/Documents/Pinpin_Clonapin/HeightPinpinClonapin/data_DRYAD/GenomicData_5165SNPs_523clones.csv", row.names=1)




# 1/ Identifying the global PEAs ####
# =================================="

# We load the piMASS outputs:
beta.snp <- read_delim("data_DRYAD/height_all_sites_res.mcmc.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
beta.snp <- as.data.frame(beta.snp)
row.names(beta.snp) <- beta.snp$rs
beta.snp$rs <- NULL

# Checking that the genotypes names in the two files are the same and are ranked in the same order.
compare::compare(row.names(beta.snp),row.names(geno))

# We select the 350 SNPs that show the strongest association with height
beta.snp.order <- beta.snp[order(abs(beta.snp$betarb),decreasing=T),]
beta.snp.order350 <- beta.snp.order[1:350,]

DF <- data.frame(gPEAs=row.names(beta.snp.order350),rPEAs_FrenchAtlantic=NA,rPEAs_IberianAtlantic=NA,rPEAs_Med=NA)





# 2/ Identifying the regional PEAs ####
# ==============================="

# 2.a/ French Atlantic region ####
# ================================"

# We load the piMASS outputs:
beta.snp <- read_delim("data_DRYAD/height_french_atlantic_res.mcmc.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
beta.snp <- as.data.frame(beta.snp)
row.names(beta.snp) <- beta.snp$rs
beta.snp$rs <- NULL

# We select the 350 SNPs that show the strongest association with height
beta.snp.order <- beta.snp[order(abs(beta.snp$betarb),decreasing=T),]
beta.snp.order350 <- beta.snp.order[1:350,]

DF$rPEAs_FrenchAtlantic <- row.names(beta.snp.order350)


# 2.b/ Iberian Atlantic region ####
# ================================="

beta.snp <- read_delim("data_DRYAD/height_iberian_atlantic_res.mcmc.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
beta.snp <- as.data.frame(beta.snp)
row.names(beta.snp) <- beta.snp$rs
beta.snp$rs <- NULL

beta.snp.order <- beta.snp[order(abs(beta.snp$betarb),decreasing=T),]
beta.snp.order350 <- beta.snp.order[1:350,]

DF$rPEAs_IberianAtlantic <- row.names(beta.snp.order350)


# 2.c/ Mediterranean region ####
# =============================="

beta.snp <- read_delim("data_DRYAD/height_mediterranean_res.mcmc.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
beta.snp <- as.data.frame(beta.snp)
row.names(beta.snp) <- beta.snp$rs
beta.snp$rs <- NULL

beta.snp.order <- beta.snp[order(abs(beta.snp$betarb),decreasing=T),]
beta.snp.order350 <- beta.snp.order[1:350,]

DF$rPEAs_Med <- row.names(beta.snp.order350)

write.csv(DF, file= "data/PEAsNames.csv", row.names = F)
saveRDS(DF,file="data/PEAsNames.RDS")


