##################################################################################################################
#                                                                                                                #
#                      Calculating gene pool-specific genomic relationship matrices                              #
#                                                                                                                #
#                                   Juliette Archambeau                                                          #
#                                       02/03/2022                                                               #
#                                                                                                                #
##################################################################################################################

# In this script, we calculate the gene-pool specific genomic relationship matrices (GRM) following the steps below:
# 1/ keeping only alleles that are not associated with height (i.e. the 'neutral' alleles)
# 2/ calculate a global GRM
# 3/ calculate gene-pool specific GRM

# More details can be found here:
#  - for steps 1 and 2: https://github.com/JulietteArchambeau/HeightPinpinClonapin/blob/master/reports/ExplanatoryVariables/GenomicRelationshipMatrices.Rmd
#  - for step 3: https://github.com/JulietteArchambeau/HeightPinpinClonapin/blob/master/reports/ExplanatoryVariables/GenePoolSpecificKinshipMatrix.Rmd


# Packages:
library(readr)
library(AGHmatrix)
library(gplots)
library(matrixcalc)
library(Matrix)
library(tidyverse)
library(bdsmatrix)



# Load the genomic data:
geno <- read.csv("~/Documents/Pinpin_Clonapin/HeightPinpinClonapin/data_DRYAD/GenomicData_5165SNPs_523clones.csv", row.names=1)



# 1/ Splitting 'neutral' and height-associated SNPs
# =================================================


# We use the outputs from a previous study: de Miguel et al. 2022.
# Here the link to the study: https://onlinelibrary.wiley.com/doi/full/10.1111/mec.16367?casa_token=1nNTc88Iy40AAAAA%3ALd4EOK5ehk_cEHIkw5A9l8nk0NPzUzlYPX8eAjVCikIjHP0WJ1kxoHJSZjMLFsZcP-8wdbNuNrlOfp1jzw
# In this study, the authors used the implemented in Bayesian variable selection regression implemented in the piMASS software to identify the SNPs associated with height
# in the five CLONAPIN common gardens. They found about 350 height associated SNPs. 
# Here we use the piMASS outputs from de Miguel et al. (2022) to remove the SNPs potentially associated with height
# and calculate the GRM only based on 'neutral' SNPs.
# piMASS = Posterior inference using Model Averaging and Subset Selection. Details here: https://stephenslab.uchicago.edu/software.html#pimass

# We load the piMASS outputs:
beta.snp <- read_delim("data/piMASSoutputs/height_all_sites_res.mcmc.txt", "\t", escape_double = FALSE, trim_ws = TRUE)



# To select about 350 SNPs considered to be associated with height,
# we set of threshold of 0.0006 for the absolute values of the Rao-Blackwellized estimates of the posterior effect size (column 'betarb')
# with this threshold, we identify 322 SNPs associated with height.

# 322 height-associated SNPs
snps.h <- beta.snp$rs[abs(beta.snp$betarb)>0.006|abs(beta.snp$betarb)==0.006]

# 4,843 SNPs not associated with height ('neutral')
snps.n <- beta.snp$rs[abs(beta.snp$betarb)<0.006] 


# We subset the genomic dataset 'geno' to obtain a dataset with only the genotypes of SNPs not associated with height.
geno <- geno[row.names(geno) %in% snps.n,] # dataset with 4843 SNPs and 523 clones



# 2/ Estimating the GRM based on all gene pools
# =============================================

# We use the function 'Gmatrix' of the package `AGHmatrix` to calculate the GRM with the VanRaden method.

# we format the genomic data for the 'Gmatrix' function:
geno[,1:dim(geno)[[2]]] <- apply(geno[ , 1:dim(geno)[[2]]], 2, function(x) as.numeric(as.character(x))) # we replace integers by numeric values
mat <- t(as.matrix(geno)) # matrix with SNPs in columns and clones in rows

# estimate the GRM:
GRM <- Gmatrix(mat,method="VanRaden", missingValue = NA, verify.posdef=T)

# to vizualize the GRM:
heatmap(GRM)
heatmap.2(as.matrix(GRM),scale="row",col=brewer.pal(11,"RdBu"),trace="none")

# The matrix has to be positive definite
# but this is not the case, 1 eigen value is lower than 0
sum(eigen(GRM)$values<0)
eigen(GRM)$values[eigen(GRM)$values<0]

# so we use of the 'nearPD' function of the package 'Matrix' to correct it:
GRM <- as.matrix(nearPD(GRM)$mat)
is.positive.definite(GRM) # now that's ok, the matrix is positive definite


# 3/ Estimating gene pool-specific GRM
# ====================================

# We follow the methodology of Muff et al. 2019
# Here the link to the paper: https://gsejournal.biomedcentral.com/track/pdf/10.1186/s12711-019-0449-7.pdf


# 3.a/ Building the Q-matrix
# ==========================


# We need the Q matrix describing the genetic population structure,
# i.e. the proportion of assignement of each clone to each gene pool

# we load the dataset containing the population structure data and keep only one row per clone
Qmat <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv") %>% as.data.frame()
Qmat <- Qmat[,c("clon",paste0("Q",rep(1:6)))]
Qmat <- unique(Qmat)

# Due to approximations, there are 18 negative values (with a value of -0.001) in the gene pool Q6
sum(Qmat[,c(paste0("Q",rep(1:6)))]<0)
filter_at(Qmat,c(paste0("Q",rep(1:6))),any_vars(. < 0))

# we set these negative values to 0.
Qmat$Q6[Qmat$Q6<0] <- 0

# clones as rownames instead of in the first column
row.names(Qmat) <- Qmat$clon
Qmat$clon <- NULL


# 3.b/ Generalized Cholesky decomposition of the GRM
# ==================================================

A.gchol <- gchol(GRM)

# Matrix T
T <- as.matrix(A.gchol)
# lower triangular matrix with diagonal elements equal to 1
# and elements below the diagonal in the respective column correspond to the expected proportion of the genome that is shared among clones.

# Transpose of matrix T
Tt <- t(as.matrix(A.gchol))

# Diagonal matrix D
diag.A <- diag(A.gchol) # vector of numeric values
D <- Diagonal(x=diag.A)


# 3.c/ Gene pool-specific matrices
# ================================

diij <- rep(1,523)
Dj <- Diagonal(x=diij)

# GRM specific to the Q1 gene pool, i.e. the Northern Africa gene pool
D1 <- Diagonal(x=Qmat[,1])
T1 <- T %*% D1
A1 <- T1 %*%  Dj %*% t(T1)
is.symmetric.matrix(as.matrix(A1))  # The matrix is symetric
is.positive.definite(as.matrix(A1)) # The matrix is not positive definite
A1 <- nearPD(A1)                    # we approximate the matrix to the nearest positive definite matrix
write.csv(as.matrix(A1$mat), file= paste0("data_DRYAD/GRM_A1.csv"))


# We do the same for the other gene pools

# Q2 gene pool, i.e. the Corsican gene pool
D2 <- Diagonal(x=Qmat[,2])
T2 <- T %*% D2
A2 <- T2 %*%  Dj %*% t(T2)
is.symmetric.matrix(as.matrix(A2))
is.positive.definite(as.matrix(A2))
A2 <- nearPD(A2)
write.csv(as.matrix(A2$mat), file= paste0("data_DRYAD/GRM_A2.csv"))

# Q3 gene pool, i.e. the Central Spain gene pool
D3 <- Diagonal(x=Qmat[,3])
T3 <- T %*% D3
A3 <- T3 %*%  Dj %*% t(T3)
is.symmetric.matrix(as.matrix(A3))
is.positive.definite(as.matrix(A3))
A3 <- nearPD(A3)
write.csv(as.matrix(A3$mat), file= paste0("data_DRYAD/GRM_A3.csv"))


# Q4 gene pool, i.e. the Atlantic French gene pool
D4 <- Diagonal(x=Qmat[,4])
T4 <- T %*% D4
A4 <- T4 %*%  Dj %*% t(T4)
is.symmetric.matrix(as.matrix(A4))
is.positive.definite(as.matrix(A4))
A4 <- nearPD(A4)
write.csv(as.matrix(A4$mat), file= paste0("data_DRYAD/GRM_A4.csv"))


# Q5 gene pool, i.e. the Atlantic Iberian gene pool
D5 <- Diagonal(x=Qmat[,5])
T5 <- T %*% D5
A5 <- T5 %*%  Dj %*% t(T5)
is.symmetric.matrix(as.matrix(A5))
is.positive.definite(as.matrix(A5))
A5 <- nearPD(A5)
write.csv(as.matrix(A5$mat), file= paste0("data_DRYAD/GRM_A5.csv"))


# Q6 gene pool, i.e. the south-eastern Spain gene pool
D6 <- Diagonal(x=Qmat[,6])
T6 <- T %*% D6
A6 <- T6 %*%  Dj %*% t(T6)
is.symmetric.matrix(as.matrix(A6))
is.positive.definite(as.matrix(A6))
A6 <- nearPD(A6)
write.csv(as.matrix(A6$mat), file= paste0("data_DRYAD/GRM_A6.csv"))
