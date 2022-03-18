########################################################################################################################"
#                                                                                                                      #
#                                     Qst-Fst analysis                                                                 #
#                                                                                                                      #
#                                   Juliette Archambeau                                                                #
#                                       18/03/2022                                                                     #
#                                                                                                                      #
########################################################################################################################"

# To determine whether height growth shows footprints of adaptive differentiation, we performed a Qst - Fst analysis.
# We used the global Fst estimate calculated in de Miguel et al. (2022) on the same data as our study (i.e. the 5,165 SNPs from
# the Illumina Infinium SNP array). de Miguel et al. (2022) used the diveRsity R package and 1,000 bootstrap iterations across
# loci to estimate the 95% confidence interval of the global Fst . They obtained a Fst of 0.112 (95% confidence interval: 0.090-0.141).

# To calculate the Qst , we used the following formula from Spitze (1993): 
# Qst = Var(provenance) / ( Var(provenance) + 2 * Var(clones) ) 
# where Var(provenance) is the variance among provenances, 
# and Var(clones) is the variance among clones (.i.e. genotypes) within provenances.

# Load the M1 model fitted on the P1 partition:
mod <- readRDS(file="outputs/models/P1/MOD1.rds")

# Extract the standard deviation of the provenances and clones
sd_prov <- as.array(mod, pars = c("sd_prov__Intercept"), fixed = TRUE) # provenances
sd_clon <- as.array(mod, pars = c("sd_prov:clon__Intercept"), fixed = TRUE) # clones

# Get the variances
sigma_prov <- sd_prov*sd_prov
sigma_clon <- sd_clon*sd_clon

# Calculate the Qst
qst <- sigma_prov / (sigma_prov + 2 * sigma_clon)

# Extract the 95% confidence interval
prob=0.95
probs <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
quantile(qst, probs = probs)

# Extract the median
median(qst)
