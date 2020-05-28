########################################################################################################"
##################                                                               #######################"
##################              Models fitted on the P2 partition                #######################"
##################                  with the package "brms"                      #######################"
##################                                                               #######################"
########################################################################################################"

# P2 partition:
# Randomly sampled provenances
# 28 provenances in the train dataset


# LIBRARIES:
library(dplyr)
library(brms)
options(mc.cores = parallel::detectCores())


# >> DATA ####
data <- readRDS(file="data/TrainP2.RDS")

# AGE
#
data$age.sc <- (data$age - mean(data$age)) / sd(data$age)


## Population structure
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


## Site climatic similarity
# Load the variance-covariance matrix of env. variables
envmat_site <- readRDS(file="data/EnvMat/site_sample1_50percent_varmat.rds")
envmat <- readRDS(file="data/EnvMat/site_sample1_50percent_varmat.rds")
row.names(envmat_site)
# Creating a column for the random effect of env. covariance
data$site_age <- paste0(data$site,data$age)


# Provenance climatic similarity
# Load the variance-covariance matrix of env. variables
envmat_prov <- readRDS(file="data/EnvMat/prov_sample4_28provs_varmat.rds")
row.names(envmat_prov)
# Creating a column for the random effect of env. covariance between provenances
data$prov_clim <- data$prov


# Count of PEAS
snp.counts <- readRDS(file="data/CountPEAs.RDS")
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
data$ppet_max_1y_site.sc <- (data$ppet_max_1y_site - mean(data$ppet_max_1y_site)) / sd(data$ppet_max_1y_site)

# Provenance climatic variables
data$bio14_prov.sc <- (data$bio14_prov - mean(data$bio14_prov)) / sd(data$bio14_prov)
data$bio5_prov.sc <- (data$bio5_prov - mean(data$bio5_prov)) / sd(data$bio5_prov)


#****************************************************************************************************************************** ####
#           MODEL 0                                                                                                             ####

mod0 <- brm(log(height) ~  age.sc + I(age.sc^2) + (1|site/block), 
            data = data, family = "gaussian",
            
            prior = c(prior(normal(0, 1), "b"),
                      prior(normal(0, 5), "Intercept")),
            
            control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=2000)

saveRDS(mod0, file="outputs/models/P2/MOD0.rds")

# Warning messages:
#   1: There were 1 divergent transitions after warmup. Increasing adapt_delta above 0.999 may help. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
# 2: Examine the pairs() plot to diagnose sampling problems




#****************************************************************************************************************************** ####
#           MODEL 1                                                                                                             ####

mod1 <- brm(log(height) ~  age.sc + I(age.sc^2) + (1|prov/clon) + (1|site/block), 
            data = data, family = "gaussian",
            
            prior = c(prior(normal(0, 1), "b"),
                      prior(normal(0, 5), "Intercept")),
            
            control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=2000)

saveRDS(mod1, file="outputs/models/P2/MOD1.rds")

# No warnings



#****************************************************************************************************************************** ####
#           MODEL 2                                                                                                             ####

mod2 <- brm(log(height) ~  age.sc + I(age.sc^2) + (1|prov/clon) + (1|site/block) + (1|prov:site), 
            data = data, family = "gaussian", 
            
            prior = c(prior(normal(0, 1), "b"),
                      prior(normal(0, 5), "Intercept")),
            
            control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=2000)

saveRDS(mod2, file="outputs/models/P2/MOD2.rds")

# Warning messages:
# 1: There were 1 divergent transitions after warmup. Increasing adapt_delta above 0.999 may help. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
# 2: Examine the pairs() plot to diagnose sampling problems




#****************************************************************************************************************************** ####
#           MODEL 7                                                                                                             ####

mod7 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
               (1|mm(Q1,Q2,Q3,Q4,Q5,Q6, weights = cbind(prop_Q1,prop_Q2,prop_Q3,prop_Q4,prop_Q5,prop_Q6))) +
               (bio5_prov.sc + bio14_prov.sc + gPEA.sc|site) + (1|block), 
             
             prior = c(prior(normal(0, 1), "b"),
                       prior(normal(0, 5), "Intercept")),
             
             data = data, family = "gaussian",
             control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)
saveRDS(mod7, file="outputs/models/P2/MOD7.rds")



#****************************************************************************************************************************** ####
#           MODEL 8                                                                                                             ####

mod8 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
               (1|mm(Q1,Q2,Q3,Q4,Q5,Q6, weights = cbind(prop_Q1,prop_Q2,prop_Q3,prop_Q4,prop_Q5,prop_Q6))) +
               (bio5_prov.sc + bio14_prov.sc + rPEA.sc|site) + (1|block), 
             
             prior = c(prior(normal(0, 1), "b"),
                       prior(normal(0, 5), "Intercept")),
             
             data = data, family = "gaussian",
             control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)
saveRDS(mod8, file="outputs/models/P2/MOD8.rds")

# No warnings!


#****************************************************************************************************************************** ####
#           MODEL 9                                                                                                             ####

mod9 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
              (1|mm(Q1,Q2,Q3,Q4,Q5,Q6, weights = cbind(prop_Q1,prop_Q2,prop_Q3,prop_Q4,prop_Q5,prop_Q6))) +
              (1|site/block), 
            
            prior = c(prior(normal(0, 1), "b"),
                      prior(normal(0, 5), "Intercept")),
            
            data = data, family = "gaussian",
            control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)
saveRDS(mod9, file="outputs/models/P2/MOD9.rds")

# Warning messages:
#   1: There were 1 divergent transitions after warmup. Increasing adapt_delta above 0.999 may help. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
# 2: Examine the pairs() plot to diagnose sampling problems





#****************************************************************************************************************************** ####
#          MODEL 10                                                                                                             ####

mod10 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
            (bio5_prov.sc + bio14_prov.sc|site) + (1|block), 
            
            prior = c(prior(normal(0, 1), "b"),
                      prior(normal(0, 5), "Intercept")),
            
            data = data, family = "gaussian",
            control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)
saveRDS(mod10, file="outputs/models/P2/MOD10.rds")

# Warning messages:
#   1: There were 1 divergent transitions after warmup. Increasing adapt_delta above 0.999 may help. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
# 2: There were 38 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 14. See
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
saveRDS(mod11, file="outputs/models/P2/MOD11.rds")


# Warning messages:
#   1: There were 1 divergent transitions after warmup. Increasing adapt_delta above 0.999 may help. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
# 2: There were 1151 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 14. See
# http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 3: Examine the pairs() plot to diagnose sampling problems


#****************************************************************************************************************************** ####
#          MODEL 12                                                                                                             ####

mod12 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
              (rPEA.sc|site) + (1|block), 
            
            prior = c(prior(normal(0, 1), "b"),
                      prior(normal(0, 5), "Intercept")),
            
            data = data, family = "gaussian",
            control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)
saveRDS(mod12, file="outputs/models/P2/MOD12.rds")


# Warning messages:
#   1: There were 3 divergent transitions after warmup. Increasing adapt_delta above 0.999 may help. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
# 2: There were 17 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 14. See
# http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 3: Examine the pairs() plot to diagnose sampling problems