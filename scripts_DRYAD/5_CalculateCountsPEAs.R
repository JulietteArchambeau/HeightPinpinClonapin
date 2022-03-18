########################################################################################################################"
#                                                                                                                      #
#                  Calculating counts of height-associated positive-effect alleles (PEAs)                              #
#                                                                                                                      #
#                                   Juliette Archambeau                                                                #
#                                       18/03/2022                                                                     #
#                                                                                                                      #
########################################################################################################################"

# In this script, we calculate counts of global and regional height-associated positive-effect alleles, rPEAs and gPEAs respectively (sections 1 and 2). 
# In section 3, we build Venn diagrams to visualize the shared proportion of gPEAs and rPEAs among regions 
# => figures in the section 2.2 of the Supplementary Information in the Am Nat manuscript

# We use the outputs of four GWAS performed in de Miguel et al. (2022) based on the same phenotypes and the same SNP marker set as in the present study
# In de Miguel et al. (2022), the authors used the Bayesian variable selection regression (BVSR) methodology
# which is implemented in the piMASS software (Guan and Stephens 2011) and corrects for population structure.
# In the four GWAS, the authors used height BLUPs accounting for site and block effects.

# First, a global GWAS was performed to identify SNPs that have an association with height at range-wide geographical scales, 
# thus using the combined phenotypic data from the five common gardens.
# and the ouputs of this global GWAS are in the file 'height_all_sites_res.mcmc.txt'.

# Second, three regional GWAS were performed to identify SNPs that have a local association with height 
# in a particular geographical region (i.e. in a particular environment), 
# thus using separately data from:
# - the Iberian Atlantic common gardens (Asturias and Portugal) => outputs in the file 'height_iberian_atlantic_res.mcmc.txt'
# - the French Atlantic common garden (Bordeaux) => outputs in the file 'height_french_atlantic_res.mcmc.txt'
# - the Mediterranean common gardens (Madrid and Cáceres) => outputs in the file 'height_mediterranean_res.mcmc.txt'

# In this script, for each of the four GWAS, we select the 350 SNPs (∼7% top associations) with the highest absolute Rao-Blackwellized estimates of the posterior effect size,
# corresponding approximately to the estimated number of SNPs with non-zero effects on height in de Miguel et al. (2022).
# These SNPs are then used to compute the counts of global and regional positive-effect alleles (gPEAs and rPEAs) for each genotype.

# More details can be found here: https://github.com/JulietteArchambeau/HeightPinpinClonapin/blob/master/reports/ExplanatoryVariables/CalculatingGPEAandRPEA.Rmd


# Packages:
library(ggvenn)  # CRAN v0.1.9
library(cowplot) # CRAN v1.0.0
library(readr)   # CRAN v1.3.1

# Load the genomic data:
geno <- read.csv("~/Documents/Pinpin_Clonapin/HeightPinpinClonapin/data_DRYAD/GenomicData_5165SNPs_523clones.csv", row.names=1)


# 1/ Counts of global PEAs ####
# ============================="

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

# we separate the SNPs that show a positive or negative association with height
snps.h.pos.350 <- row.names(beta.snp.order350[beta.snp.order350$betarb>0,]) # 186 SNPs with a positive association
snps.h.neg.350 <- row.names(beta.snp.order350[beta.snp.order350$betarb<0,]) # 164 SNPs with a negative association

# We create a list to compare PEA overlapp and direction of their effect in section 3.
pea.list <- list()
pea.list$gPEAs <- list(all=rownames(beta.snp.order350),
                       pos=snps.h.pos.350,
                       neg=snps.h.neg.350)


# In the genomic data, we subset:
geno.snp.h.pos.350 <- geno[row.names(geno) %in% snps.h.pos.350,] # SNPs positively associated wih height
geno.snp.h.neg.350 <- geno[row.names(geno) %in% snps.h.neg.350,] # SNPs negatively associated with height

# We invert 0 and 2 in the subset of alleles with negative effects, to only have alleles with a positive effect
snpnames <- row.names(geno.snp.h.neg.350)
geno.snp.h.neg.350 <- geno.snp.h.neg.350 %>%  mutate_if(is.integer, as.numeric) %>% mutate_all(list(~recode(., `0` = 2,`2`=0)))
row.names(geno.snp.h.neg.350) <- snpnames 

# we recombine positive and negative alleles
geno.snp.h.350 <-rbind(geno.snp.h.pos.350,geno.snp.h.neg.350)
names.all <- rownames(geno.snp.h.350)

# we create a df with one row per clone and in which we are going to add the counts of gPEAs and rPEAs for each clone:
DF <- data.frame(clon=names(geno.snp.h.350),count_all_350=apply(geno.snp.h.350,2,sum,na.rm=T))


# 2/ Counts of regional PEAs ####
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

# we separate the SNPs that show a positive or negative association with height
snps.h.pos.350 <- row.names(beta.snp.order350[beta.snp.order350$betarb>0,]) # 177 SNPs positively associated with height
snps.h.neg.350 <- row.names(beta.snp.order350[beta.snp.order350$betarb<0,]) # 173 SNPs negatively associated with height

# Put them in the list for looking at overlapp and direction of effect in section 3.
pea.list$rPEAs_FA <- list(all=rownames(beta.snp.order350),
                          pos=snps.h.pos.350,
                          neg=snps.h.neg.350)

# In the genomic data, we subset:
geno.snp.h.pos.350 <- geno[row.names(geno) %in% snps.h.pos.350,] # the 177 SNPs positively associated wih height
geno.snp.h.neg.350 <- geno[row.names(geno) %in% snps.h.neg.350,] # the 173 SNPs negatively associated with height

# we invert 0 and 2 in the subset of alleles with negative effects, to only have alleles with a positive effect
snpnames <- row.names(geno.snp.h.neg.350)
geno.snp.h.neg.350 <- geno.snp.h.neg.350 %>%  mutate_if(is.integer, as.numeric) %>% mutate_all(list(~recode(., `0` = 2,`2`=0)))
row.names(geno.snp.h.neg.350) <- snpnames 


# we recombine positive and negative alleles
geno.snp.h.350 <-rbind(geno.snp.h.pos.350,geno.snp.h.neg.350)
names.fr_atl <- rownames(geno.snp.h.350)

# add the counts of rPEAs to the df with the counts of gPEAs for each clone
df <- data.frame(clon=names(geno.snp.h.350),count_fratl_350=apply(geno.snp.h.350,2,sum,na.rm=T))
DF <- merge(DF,df,by="clon")


# 2.b/ Iberian Atlantic region ####
# ================================="

# Same code as in the 2.a subsection:

beta.snp <- read_delim("data_DRYAD/height_iberian_atlantic_res.mcmc.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
beta.snp <- as.data.frame(beta.snp)
row.names(beta.snp) <- beta.snp$rs
beta.snp$rs <- NULL

beta.snp.order <- beta.snp[order(abs(beta.snp$betarb),decreasing=T),]
beta.snp.order350 <- beta.snp.order[1:350,]

snps.h.pos.350 <- row.names(beta.snp.order350[beta.snp.order350$betarb>0,]) # 169 SNPs positively associated with height
snps.h.neg.350 <- row.names(beta.snp.order350[beta.snp.order350$betarb<0,]) # 181 SNPs negatively associated with height

pea.list$rPEAs_IA <- list(all=rownames(beta.snp.order350),
                          pos=snps.h.pos.350,
                          neg=snps.h.neg.350)

geno.snp.h.pos.350 <- geno[row.names(geno) %in% snps.h.pos.350,]
geno.snp.h.neg.350 <- geno[row.names(geno) %in% snps.h.neg.350,]

snpnames <- row.names(geno.snp.h.neg.350)
geno.snp.h.neg.350 <- geno.snp.h.neg.350 %>%  mutate_if(is.integer, as.numeric) %>% mutate_all(list(~recode(., `0` = 2,`2`=0)))
row.names(geno.snp.h.neg.350) <- snpnames 

geno.snp.h.350 <-rbind(geno.snp.h.pos.350,geno.snp.h.neg.350)
names.ib_atl <- rownames(geno.snp.h.350)

df <- data.frame(clon=names(geno.snp.h.350),count_ibatl_350=apply(geno.snp.h.350,2,sum,na.rm=T))

DF <- merge(DF,df,by="clon")


# 2.c/ Mediterranean region ####
# =============================="

# Same code as in the 2.a subsection:

beta.snp <- read_delim("data_DRYAD/height_mediterranean_res.mcmc.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
beta.snp <- as.data.frame(beta.snp)
row.names(beta.snp) <- beta.snp$rs
beta.snp$rs <- NULL

beta.snp.order <- beta.snp[order(abs(beta.snp$betarb),decreasing=T),]
beta.snp.order350 <- beta.snp.order[1:350,]

snps.h.pos.350 <- row.names(beta.snp.order350[beta.snp.order350$betarb>0,]) # 174 SNPs positively associated with height
snps.h.neg.350 <- row.names(beta.snp.order350[beta.snp.order350$betarb<0,]) # 176 SNPs negatively associated with height

pea.list$rPEAs_Med <- list(all=rownames(beta.snp.order350),
                           pos=snps.h.pos.350,
                           neg=snps.h.neg.350)

geno.snp.h.pos.350 <- geno[row.names(geno) %in% snps.h.pos.350,]
geno.snp.h.neg.350 <- geno[row.names(geno) %in% snps.h.neg.350,]

snpnames <- row.names(geno.snp.h.neg.350)
geno.snp.h.neg.350 <- geno.snp.h.neg.350 %>%  mutate_if(is.integer, as.numeric) %>% mutate_all(list(~recode(., `0` = 2,`2`=0)))
row.names(geno.snp.h.neg.350) <- snpnames 

geno.snp.h.350 <-rbind(geno.snp.h.pos.350,geno.snp.h.neg.350)
names.med <- rownames(geno.snp.h.350)

df <- data.frame(clon=names(geno.snp.h.350),count_med_350=apply(geno.snp.h.350,2,sum,na.rm=T))

DF <- merge(DF,df,by="clon")

# We save the dataframe with the counts of gPEAs (1 column: 'count_all_350') 
# and the counts of rPEAs (3 columns: 'count_fratl_350', 'count_ibatl_350' and 'count_med_350') 
write.csv(DF, file= "data_DRYAD/CountPEAs.csv")



# 3/ Shared proportion of gPEAs and rPEAs among regions ####
# =========================================================="

# 3.a/ Venn diagrams with only rPEAs ####
# ======================================="

pos.list <- list("Mediterranean"=pea.list$rPEAs_Med$pos,
                 "French Atlantic"=pea.list$rPEAs_FA$pos,
                 "Iberian Atlantic"=pea.list$rPEAs_IA$pos)

neg.list <- list("Mediterranean"=pea.list$rPEAs_Med$neg,
                 "French Atlantic"=pea.list$rPEAs_FA$neg,
                 "Iberian Atlantic"=pea.list$rPEAs_IA$neg)

all.list <- list("Mediterranean"=pea.list$rPEAs_Med$all,
                 "French Atlantic"=pea.list$rPEAs_FA$all,
                 "Iberian Atlantic"=pea.list$rPEAs_IA$all)


p.all <- ggvenn(all.list)
p.pos <- ggvenn(pos.list)
p.neg <- ggvenn(neg.list)

# Figure S2:
p <- plot_grid(p.all,p.pos,p.neg,
               nrow=1,
               labels=c("All SNPs","Positive-effect SNPs","Negative-effect SNPs"),
               label_size = 20)


# 3.b/ Venn diagrams with rPEAs and gPEAs ####
# ============================================"

pos.list <- list("Mediterranean"=pea.list$rPEAs_Med$pos,
                 "French Atlantic"=pea.list$rPEAs_FA$pos,
                 "Iberian Atlantic"=pea.list$rPEAs_IA$pos,
                 "Global"=pea.list$gPEAs$pos)

neg.list <- list("Mediterranean"=pea.list$rPEAs_Med$neg,
                 "French Atlantic"=pea.list$rPEAs_FA$neg,
                 "Iberian Atlantic"=pea.list$rPEAs_IA$neg,
                 "Global"=pea.list$gPEAs$neg)

all.list <- list("Mediterranean"=pea.list$rPEAs_Med$all,
                 "French Atlantic"=pea.list$rPEAs_FA$all,
                 "Iberian Atlantic"=pea.list$rPEAs_IA$all,
                 "Global"=pea.list$gPEAs$all)


p.all <- ggvenn(all.list,set_name_size=5)
p.pos <- ggvenn(pos.list,set_name_size=5)
p.neg <- ggvenn(neg.list,set_name_size=5)

# Figure S3:
p <- plot_grid(p.all,p.pos,p.neg,
               nrow=1,
               labels=c("All SNPs","Positive-effect SNPs","Negative-effect SNPs"),
               label_size = 20)

