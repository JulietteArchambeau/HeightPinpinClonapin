########################################################################################################"
##################                                                               #######################"
##################              Models fitted on the P1 partition                #######################"
##################                  with the package "brms"                      #######################"
##################                                                               #######################"
########################################################################################################"

# P1 partition:
# Randomly sampled observations
# 75% of the observations in the train dataset


# LIBRARIES:
library(dplyr)
library(brms)
options(mc.cores = parallel::detectCores())


# >> DATA ####
data <- readRDS(file="data/TrainP1.RDS")

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


# 6 weighed Genomic relationship matrices
A1 <- readRDS(file="data/GRMs/a1_grm.rds")
A2 <- readRDS(file="data/GRMs/a2_grm.rds")
A3 <- readRDS(file="data/GRMs/a3_grm.rds")
A4 <- readRDS(file="data/GRMs/a4_grm.rds")
A5 <- readRDS(file="data/GRMs/a5_grm.rds")
A6 <- readRDS(file="data/GRMs/a6_grm.rds")

# Intercepts of the 6 weighed Genomic relationship matrices
data$clon1 <- data$clon
data$clon2 <- data$clon
data$clon3 <- data$clon
data$clon4 <- data$clon
data$clon5 <- data$clon
data$clon6 <- data$clon


## Site climatic similarity
# Load the variance-covariance matrix of env. variables
envmat_site <- readRDS(file="data/EnvMat/site_sample1_75percent_varmat.rds")
row.names(envmat_site)
# Creating a column for the random effect of env. covariance
data$site_age <- paste0(data$site,data$age)


# Prov climatic similarity
# Load the variance-covariance matrix of env. variables
envmat_prov <- readRDS(file="data/EnvMat/prov_sample1_75percent_varmat.rds")
row.names(envmat_prov)
# Creating a column for the random effect of env. covariance between provenances
data$prov_clim <- data$prov


# Count of the PEAs
snp.counts <- readRDS(file="data/CountPEAs.RDS")
data <- merge(data,snp.counts,by="clon")
rm(snp.counts)

# Region-specific PEAs
data$count_all_350 <- NA 
data$count_all_350[data$site=="asturias"|data$site=="portugal"] <- data$count_ibatl_350[data$site=="asturias"|data$site=="portugal"]
data$count_all_350[data$site=="caceres"|data$site=="madrid"] <- data$count_med_350[data$site=="madrid"|data$site=="caceres"]
data$count_all_350[data$site=="bordeaux"] <- data$count_fratl_350[data$site=="bordeaux"]
data$count_all_350.sc <- scale(data$count_all_350)


# Global PEAs
data$count_all.sc <- scale(data$count_all)
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
#           MODEL 1                                                                                                             ####

mod1 <- brm(log(height) ~  age.sc + I(age.sc^2) + (1|prov/clon) + (1|site/block), 
            data = data, family = "gaussian",
            
            prior = c(prior(normal(0, 1), "b"),
              prior(normal(0, 5), "Intercept")),
            
            control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=2000)

saveRDS(mod1, file="outputs/models/P1/MOD1.rds")

# No warnings.




#****************************************************************************************************************************** ####
#           MODEL 2                                                                                                             ####

mod2 <- brm(log(height) ~  age.sc + I(age.sc^2) + (1|prov/clon) + (1|site/block) + (1|prov:site), 
                  data = data, family = "gaussian", 
                  
                  prior = c(prior(normal(0, 1), "b"),
                            prior(normal(0, 5), "Intercept")),
                  
                  control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=2000)

saveRDS(mod2, file="outputs/models/P1/MOD2.rds")

# No warnings.




#****************************************************************************************************************************** ####
#           MODEL 3                                                                                                             ####


mod3 <- brm(log(height) ~  age.sc + I(age.sc^2) + (1|clon) + (1|prov) + (1|site/block) + (1|site_age), 
            data = data, family = "gaussian", 
            
            prior = c(prior(normal(0, 1), "b"),
                      prior(normal(0, 5), "Intercept")),
            
            cov_ranef = list(site_age = envmat_site),
            
            control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=2000)
saveRDS(mod3, file="outputs/models/P1/MOD3.rds")

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

saveRDS(mod4, file="outputs/models/P1/MOD4.rds")

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



saveRDS(mod5, file="outputs/models/P1/MOD5.rds")

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
saveRDS(mod6, file="outputs/models/P1/MOD6.rds")

# No warnings.





#****************************************************************************************************************************** ####
#           MODEL 7                                                                                                             ####

mod7 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
               (1|mm(Q1,Q2,Q3,Q4,Q5,Q6, weights = cbind(prop_Q1,prop_Q2,prop_Q3,prop_Q4,prop_Q5,prop_Q6))) +
               (bio5_prov.sc + bio14_prov.sc + count_all.sc|site) + (1|block), 
             
             prior = c(prior(normal(0, 1), "b"),
                       prior(normal(0, 5), "Intercept")),
             
             data = data, family = "gaussian",
             control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)
saveRDS(mod7, file="outputs/models/P1/MOD7.rds")
# 
# Warning messages:
#   1: There were 1 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 14. See
# http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 2: Examine the pairs() plot to diagnose sampling problems





#****************************************************************************************************************************** ####
#           MODEL 8                                                                                                             ####


mod8 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
               (1|mm(Q1,Q2,Q3,Q4,Q5,Q6, weights = cbind(prop_Q1,prop_Q2,prop_Q3,prop_Q4,prop_Q5,prop_Q6))) +
               (bio5_prov.sc + bio14_prov.sc + count_all_350.sc|site) + (1|block), 
             
             prior = c(prior(normal(0, 1), "b"),
                       prior(normal(0, 5), "Intercept")),
             
             data = data, family = "gaussian",
             control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)
saveRDS(mod8, file="outputs/models/P1/MOD8.rds")

# 
# Warning messages:
#   1: There were 2 divergent transitions after warmup. Increasing adapt_delta above 0.999 may help. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
# 2: There were 1 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 14. See
# http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 3: Examine the pairs() plot to diagnose sampling problems





#****************************************************************************************************************************** ####
#           MODEL 9                                                                                                             ####

mod9 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
              (1|mm(Q1,Q2,Q3,Q4,Q5,Q6, weights = cbind(prop_Q1,prop_Q2,prop_Q3,prop_Q4,prop_Q5,prop_Q6))) +
              (1|site/block), 
            
            prior = c(prior(normal(0, 1), "b"),
                      prior(normal(0, 5), "Intercept")),
            
            data = data, family = "gaussian",
            control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)
saveRDS(mod9, file="outputs/models/P1/MOD9.rds")
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
saveRDS(mod10, file="outputs/models/P1/MOD10.rds")


# Warning messages:
#   1: There were 2 divergent transitions after warmup. Increasing adapt_delta above 0.999 may help. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
# 2: There were 9 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 14. See
# http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 3: Examine the pairs() plot to diagnose sampling problems



#****************************************************************************************************************************** ####
#          MODEL 11                                                                                                             ####

mod11 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
               (count_all.sc|site) + (1|block), 
             
             prior = c(prior(normal(0, 1), "b"),
                       prior(normal(0, 5), "Intercept")),
             
             data = data, family = "gaussian",
             control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)
saveRDS(mod11, file="outputs/models/P1/MOD11.rds")
# 
# Warning messages:
#   1: There were 1 divergent transitions after warmup. Increasing adapt_delta above 0.999 may help. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
# 2: Examine the pairs() plot to diagnose sampling problems





#****************************************************************************************************************************** ####
#          MODEL 12                                                                                                             ####

mod12 <- brm(log(height) ~  age.sc + I(age.sc^2) + 
               (count_all_350.sc|site) + (1|block), 
             
             prior = c(prior(normal(0, 1), "b"),
                       prior(normal(0, 5), "Intercept")),
             
             data = data, family = "gaussian",
             control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)
saveRDS(mod12, file="outputs/models/P1/MOD12.rds")

# Warning messages:
#   1: There were 1 divergent transitions after warmup. Increasing adapt_delta above 0.999 may help. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 
# 2: Examine the pairs() plot to diagnose sampling problems




#****************************************************************************************************************************** ####
#          MODEL 13 (model 3bis in the SuppInfo)                                                                                ####

mod13 <- brm(log(height) ~  age.sc + I(age.sc^2) + (1|site_age) +  (1|block) + (1|prov/clon), 
             
             prior = c(prior(normal(0, 1), "b"),
                       prior(normal(0, 5), "Intercept")),
             
             cov_ranef = list(site_age = envmat_site),
              
             data = data, family = "gaussian",
             control = list(adapt_delta=0.999,max_treedepth =14),chain=4,iter=3000)
saveRDS(mod13, file="outputs/models/P1/MOD13.rds")

# No warnings!