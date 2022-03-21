########################################################################################################################"
#                                                                                                                      #
#                   Fitting models on the P1 partition with the brms package                                           #
#                                                                                                                      #
#                                   Juliette Archambeau                                                                #
#                                       18/03/2022                                                                     #
#                                                                                                                      #
########################################################################################################################"

# P1 partition:
# Randomly sampled observations
# 75% of the observations in the train dataset

# All models were saved after fitting, but are not included in the DRYAD repository as they are very heavy.

# Packages:
library(dplyr) # CRAN v1.0.0
library(brms)  # CRAN v2.11.1
options(mc.cores = parallel::detectCores())
library(readr) # CRAN v1.3.1

# Load data:
data <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv") %>% 
  dplyr::filter(P1=="train")


# AGE
data$age.sc <- (data$age - mean(data$age)) / sd(data$age)


# Population structure
colnames(data)[colnames(data) %in% c(paste0("Q",rep(1:6)))] <- c(paste0("prop_Q",rep(1:6)))
head(data[,c(paste0("prop_Q",rep(1:6)))])
data$Q1 <- "Q1"
data$Q2 <- "Q2"
data$Q3 <- "Q3"
data$Q4 <- "Q4"
data$Q5 <- "Q5"
data$Q6 <- "Q6"
sum(data[,c(paste0("prop_Q",rep(1:6)))]<0)
filter_at(data,c(paste0("prop_Q",rep(1:6))),any_vars(. < 0))
data$prop_Q6[data$prop_Q6<0] <- 0


# 6 gene pool-specific genomic relationship matrices
A1 <- read.csv("data_DRYAD/GRM_A1.csv", row.names=1) %>% as.matrix()
A2 <- read.csv("data_DRYAD/GRM_A2.csv", row.names=1) %>% as.matrix()
A3 <- read.csv("data_DRYAD/GRM_A3.csv", row.names=1) %>% as.matrix()
A4 <- read.csv("data_DRYAD/GRM_A4.csv", row.names=1) %>% as.matrix()
A5 <- read.csv("data_DRYAD/GRM_A5.csv", row.names=1) %>% as.matrix()
A6 <- read.csv("data_DRYAD/GRM_A6.csv", row.names=1) %>% as.matrix()

# Intercepts of the 6 GRMs
data$clon1 <- data$clon
data$clon2 <- data$clon
data$clon3 <- data$clon
data$clon4 <- data$clon
data$clon5 <- data$clon
data$clon6 <- data$clon


## Site climatic similarity
# Load the variance-covariance matrix of env. variables
envmat_site <- read.csv(file="data_DRYAD/VarCovMatSites.csv", row.names=1) %>% as.matrix()
               # Creating a column for the random effect of env. covariance
data$site_age <- paste0(data$site,data$age)


# Provenance climatic similarity
# Load the variance-covariance matrix of env. variables
envmat_prov <- read.csv(file="data_DRYAD/VarCovMatProvenancesP1.csv", row.names=1) %>% as.matrix()
# Creating a column for the random effect of env. covariance between provenances
data$prov_clim <- data$prov


# Count of the PEAs
snp.counts <- read.csv(file="data_DRYAD/CountPEAs.csv", row.names = 1)
data <- merge(data,snp.counts,by="clon")
rm(snp.counts)

# Region-specific PEAs
data$rPEA <- NA 
data$rPEA[data$site=="asturias"|data$site=="portugal"] <- data$count_ibatl_350[data$site=="asturias"|data$site=="portugal"]
data$rPEA[data$site=="caceres"|data$site=="madrid"] <- data$count_med_350[data$site=="madrid"|data$site=="caceres"]
data$rPEA[data$site=="bordeaux"] <- data$count_fratl_350[data$site=="bordeaux"]
data$rPEA.sc <- scale(data$rPEA)


# Global PEAs
data$gPEA <- data$count_all_350
data$gPEA.sc <- scale(data$gPEA)
names(data)



# Site climatic variables
data$pre_mean_1y_site.sc <- (data$pre_mean_1y_site - mean(data$pre_mean_1y_site)) / sd(data$pre_mean_1y_site)
data$pre_summer_min_site.sc <- (data$pre_summer_min_site - mean(data$pre_summer_min_site)) / sd(data$pre_summer_min_site)
data$tmn_min_1y_site.sc <- (data$tmn_min_1y_site - mean(data$tmn_min_1y_site)) / sd(data$tmn_min_1y_site)
data$tmx_max_1y_site.sc <- (data$tmx_max_1y_site - mean(data$tmx_max_1y_site)) / sd(data$tmx_max_1y_site)
data$tmx_mean_1y_site.sc <- (data$tmx_mean_1y_site - mean(data$tmx_mean_1y_site)) / sd(data$tmx_mean_1y_site)


# Provenance climatic variables
data$bio14_prov.sc <- (data$bio14_prov - mean(data$bio14_prov)) / sd(data$bio14_prov)
data$bio5_prov.sc <- (data$bio5_prov - mean(data$bio5_prov)) / sd(data$bio5_prov)




#****************************************************************************************************************************** ####
#           MODEL 1                                                                                                             ####

mod1 <- brm(log(height) ~  age.sc + I(age.sc^2) + (1|prov/clon) + (1|site/block), 
            data = data, family = "gaussian",
            
            prior = c(prior(normal(0, 1), "b"),
                      prior(normal(0, 5), "Intercept")),
            
            control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=2000)


# No warnings.




#****************************************************************************************************************************** ####
#           MODEL 2                                                                                                             ####

mod2 <- brm(log(height) ~  age.sc + I(age.sc^2) + (1|prov/clon) + (1|site/block) + (1|prov:site), 
            data = data, family = "gaussian", 
            
            prior = c(prior(normal(0, 1), "b"),
                      prior(normal(0, 5), "Intercept")),
            
            control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=2000)


# No warnings.




#****************************************************************************************************************************** ####
#           MODEL 3                                                                                                             ####


mod3 <- brm(log(height) ~  age.sc + I(age.sc^2) + (1|clon) + (1|prov) + (1|site/block) + (1|site_age), 
            data = data, family = "gaussian", 
            
            prior = c(prior(normal(0, 1), "b"),
                      prior(normal(0, 5), "Intercept")),
            
            cov_ranef = list(site_age = envmat_site),
            
            control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=2000)

# 
# Warning message:
#   Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#tail-ess 






#****************************************************************************************************************************** ####
#           MODEL 4                                                                                                             ####

mod4 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
              (1|mm(Q1,Q2,Q3,Q4,Q5,Q6, weights = cbind(prop_Q1,prop_Q2,prop_Q3,prop_Q4,prop_Q5,prop_Q6))) +
              (1|prov/clon) + (1|site/block) + (1|site_age),
            
            prior = c(prior(normal(0, 1), "b"),
                      prior(normal(0, 5), "Intercept")),
            
            cov_ranef = list(site_age = envmat_site),
            
            data = data, family = "gaussian",
            control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=2500)


# Warning messages:
#   1: There were 6 divergent transitions after warmup. Increasing adapt_delta above 0.999 may help. See
#   http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
#   2: Examine the pairs() plot to diagnose sampling problems





#****************************************************************************************************************************** ####
#           MODEL 5                                                                                                             ####

mod5 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
              (1|clon1) + (1|clon2) + (1|clon3) + (1|clon4) + (1|clon5) + (1|clon6) +
              (1|mm(Q1,Q2,Q3,Q4,Q5,Q6, weights = cbind(prop_Q1,prop_Q2,prop_Q3,prop_Q4,prop_Q5,prop_Q6))) +
              (1|site/block) + (1|site_age) + (1|prov),
            
            prior = c(prior(normal(0, 1), "b"),
                      prior(normal(0, 5), "Intercept")),
            
            cov_ranef = list(site_age = envmat_site,
                             clon1=A1,clon2=A2,clon3=A3,clon4=A4,clon5=A5,clon6=A6),
            
            data = data, family = "gaussian",
            control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)

# No warnings.





#****************************************************************************************************************************** ####
#           MODEL 6                                                                                                             ####

mod6 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
              (1|mm(Q1,Q2,Q3,Q4,Q5,Q6, weights = cbind(prop_Q1,prop_Q2,prop_Q3,prop_Q4,prop_Q5,prop_Q6))) +
              (1|prov/clon) + (1|site/block)  + (1|site_age) + (1|prov_clim), 
            
            prior = c(prior(normal(0, 1), "b"),
                      prior(normal(0, 5), "Intercept")),
            
            cov_ranef = list(site_age = envmat_site,prov_clim = envmat_prov),
            
            data = data, family = "gaussian",
            control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)

# No warnings.





#****************************************************************************************************************************** ####
#           MODEL 7                                                                                                             ####

mod7 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
              (1|mm(Q1,Q2,Q3,Q4,Q5,Q6, weights = cbind(prop_Q1,prop_Q2,prop_Q3,prop_Q4,prop_Q5,prop_Q6))) +
              (bio5_prov.sc + bio14_prov.sc + gPEA.sc|site) + (1|block), 
            
            prior = c(prior(normal(0, 1), "b"),
                      prior(normal(0, 5), "Intercept")),
            
            data = data, family = "gaussian",
            control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)


# No warnings!

#****************************************************************************************************************************** ####
#           MODEL 8                                                                                                             ####


mod8 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
              (1|mm(Q1,Q2,Q3,Q4,Q5,Q6, weights = cbind(prop_Q1,prop_Q2,prop_Q3,prop_Q4,prop_Q5,prop_Q6))) +
              (bio5_prov.sc + bio14_prov.sc + rPEA.sc|site) + (1|block), 
            
            prior = c(prior(normal(0, 1), "b"),
                      prior(normal(0, 5), "Intercept")),
            
            data = data, family = "gaussian",
            control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)

# Warning messages:
#   1: There were 37 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 14. See
# http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 2: Examine the pairs() plot to diagnose sampling problems
# 



#****************************************************************************************************************************** ####
#           MODEL 9                                                                                                             ####

mod9 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
              (1|mm(Q1,Q2,Q3,Q4,Q5,Q6, weights = cbind(prop_Q1,prop_Q2,prop_Q3,prop_Q4,prop_Q5,prop_Q6))) +
              (1|site/block), 
            
            prior = c(prior(normal(0, 1), "b"),
                      prior(normal(0, 5), "Intercept")),
            
            data = data, family = "gaussian",
            control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)

# 
# Warning messages:
#   1: There were 2 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 14. See
# http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 2: Examine the pairs() plot to diagnose sampling problems





#****************************************************************************************************************************** ####
#          MODEL 10                                                                                                             ####

mod10 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
               (bio5_prov.sc + bio14_prov.sc|site) + (1|block), 
             
             prior = c(prior(normal(0, 1), "b"),
                       prior(normal(0, 5), "Intercept")),
             
             data = data, family = "gaussian",
             control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)


# Warning messages:
#   1: There were 2 divergent transitions after warmup. Increasing adapt_delta above 0.999 may help. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
# 2: There were 9 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 14. See
# http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 3: Examine the pairs() plot to diagnose sampling problems



#****************************************************************************************************************************** ####
#          MODEL 11                                                                                                             ####

mod11 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
               (gPEA.sc|site) + (1|block), 
             
             prior = c(prior(normal(0, 1), "b"),
                       prior(normal(0, 5), "Intercept")),
             
             data = data, family = "gaussian",
             control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)

# Warning messages:
#   1: There were 2 divergent transitions after warmup. Increasing adapt_delta above 0.999 may help. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
# 2: There were 27 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 14. See
# http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 3: Examine the pairs() plot to diagnose sampling problems
# 
# 4: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#bulk-ess 



#****************************************************************************************************************************** ####
#          MODEL 12                                                                                                             ####

mod12 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
               (rPEA.sc|site) + (1|block), 
             
             prior = c(prior(normal(0, 1), "b"),
                       prior(normal(0, 5), "Intercept")),
             
             data = data, family = "gaussian",
             control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)

# Warning messages:
#   1: There were 2 divergent transitions after warmup. Increasing adapt_delta above 0.999 may help. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
# 2: There were 36 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 14. See
# http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 3: Examine the pairs() plot to diagnose sampling problems



#****************************************************************************************************************************** ####
#          MODEL 13 (model 3bis in the Supplementary Information, section 6.1.2)                                                ####

mod13 <- brm(log(height) ~  age.sc + I(age.sc^2) + (1|site_age) +  (1|block) + (1|prov/clon), 
             
             prior = c(prior(normal(0, 1), "b"),
                       prior(normal(0, 5), "Intercept")),
             
             cov_ranef = list(site_age = envmat_site),
             
             data = data, family = "gaussian",
             control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)

# No warnings!
