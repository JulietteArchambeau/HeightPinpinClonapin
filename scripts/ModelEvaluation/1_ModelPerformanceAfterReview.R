###############################################################################################"
##################                                                      #######################"
##################                  Model performance                   #######################"
##################                    After review                      #######################"
##################                                                      #######################"
###############################################################################################"


library(loo)
library(dplyr)
library(tidyverse)
library(stringr)
library(xtable)
library(brms)
library(gdata)


# Which baseline model? 
baseline <- "MOD2"


# Data partition used to fit the models:
part <- "P1" # choose between P1, P2 and P3

# Path
path= paste0("outputs/models/",part,"/")


# Models in a list

myFiles <- list.files(path=path,pattern=".rds")
models <- list()
for (i in 1:length(myFiles)){
  models[[i]] <- readRDS(file=paste0(path,myFiles[i]))
}
names(models) <- str_sub(myFiles,0,-5)



# Traning and test datasets:
if (grepl("P1",path)==TRUE){
  test <- readRDS(file="data/TestP1prepared.RDS")
  train <- readRDS(file="data/TrainP1.RDS")
} else if (grepl("P2",path)==TRUE){
  test <- readRDS(file="data/TestP2prepared.RDS")
  train <- readRDS(file="data/TrainP2.RDS")
} else if (grepl("P3",path)==TRUE){
  test <- readRDS(file="data/TestP3prepared.RDS")
  train <- readRDS(file="data/TrainP3.RDS")
}


# GLOBAL TABLES (all sites merged) ####
# ====================================="

# LOO    ####
# =========="

# Load loo table
myloos <- list.files(path=paste0("outputs/loos/",part,"/"),pattern=".rds") 
myloos <- str_sub(myloos,5,-5)
loos <- list()

for (i in 1:length(myloos)){
  loos[[i]] <- readRDS(file=paste0("outputs/loos/",part,"/","loo_",myloos[i],".rds"))
}
names(loos) <- myloos
df <- as.data.frame(loo::loo_compare(loos))
df






# R2 - Proportion of variance explained in the test and train datasets  #####
# =========================================================================="



# R2 conditional on age
########################"

# Variance explained by age in M1
  # In the training dataset:
var_age_baseline_train  <- pp_expect(models[[baseline]], transform = TRUE,re.form=NA) %>% 
  apply(1,var) %>% 
  mean

  # In the test dataset
var_age_baseline_test  <- pp_expect(models[[baseline]], transform = TRUE,newdata=test,allow_new_levels=T,re.form=NA) %>% 
  apply(1,var) %>% 
  mean

# Functions to calculate R2|age in the training and test datasets, respectively:
R2condAge_SameData <- function(fit,re.form) {
  y <- rstanarm::get_y(fit)
  var_y <- var(y)
  ypred <- pp_expect(fit, transform = TRUE,re.form=re.form)
  var_ypred <- apply(ypred, 1, var)
  (var_ypred - var_age_baseline_train)/ (var_y - var_age_baseline_train)
}

R2condAge_OutofSample <- function(fit,re.form) {
  y <- log(test$height)
  var_y <- var(y)
  ypred <- pp_expect(fit, transform = TRUE,re.form=re.form,newdata = test,allow_new_levels=T)
  var_ypred <- apply(ypred, 1, var)
  (var_ypred - var_age_baseline_test) / (var_y - var_age_baseline_test)
}

# Calculate R2|age in the training and test datasets, respectively:
for (i in 1:length(models)){
  R2condAge  <- R2condAge_SameData(models[[i]],re.form = NULL)  %>% 
    posterior_summary(probs=c(0.025,0.975)) %>% 
    as_tibble() 
  df[names(models[i]),"R2condAge_Train"] <- paste0(round(R2condAge$Estimate, 3)," [",
                                                   round(R2condAge$Q2.5, 3),"-",
                                                   round(R2condAge$Q97.5, 3),"]")
  
  R2condAge  <- R2condAge_OutofSample(models[[i]],re.form = NULL)  %>% 
    posterior_summary(probs=c(0.025,0.975)) %>% 
    as_tibble() 
  df[names(models[i]),"R2condAge_Test"] <- paste0(round(R2condAge$Estimate, 3)," [",
                                                  round(R2condAge$Q2.5, 3),"-",
                                                  round(R2condAge$Q97.5, 3),"]")
  
  }



# R2 in sample, i.e. proportion of variance explained on the same data (training dataset)
#########################################################################################"

# 2 values are given: the mean of which is used as the measure of central tendency and the standard deviation as the measure of variability
R2_over_VarTot_SameData <- function(fit,re.form) {
  y <- rstanarm::get_y(fit)
  var_y <- var(y)
  ypred <- pp_expect(fit, transform = TRUE,re.form=re.form)
  var_ypred <- apply(ypred, 1, var)
  var_ypred / var_y
}

for (i in 1:length(models)){
  R2all  <- R2_over_VarTot_SameData(models[[i]],re.form = NULL)  %>% 
    posterior_summary(probs=c(0.025,0.975),robust=T) %>% 
    as_tibble() 
  R2fix  <- R2_over_VarTot_SameData(models[[i]],re.form=NA)  %>% 
    posterior_summary(probs=c(0.025,0.975),robust=T) %>% 
    as_tibble() 
  df[names(models[i]),"R2all_Train"] <- paste0(round(R2all$Estimate, 3)," [",
                                               round(R2all$Q2.5, 3),"-",
                                               round(R2all$Q97.5, 3),"]") #  R2 of random + fixed effects
  df[names(models[i]),"R2fix_Train"] <- paste0(round(R2fix$Estimate, 3)," [",
                                               round(R2fix$Q2.5, 3),"-",
                                               round(R2fix$Q97.5, 3),"]")  # R2 of fixed effects
}


# R2 out-of-sample, i.e. the proportion of variance explained on new observations or new provenances (test data set)
####################################################################################################################"

R2_over_VarTot_OutofSample <- function(fit,re.form) {
  y <- log(test$height)
  var_y <- var(y)
  ypred <- pp_expect(fit, transform = TRUE,re.form=re.form,newdata = test,allow_new_levels=T)
  var_ypred <- apply(ypred, 1, var)
  var_ypred / var_y
}

for (i in 1:length(models)){
  R2all  <- R2_over_VarTot_OutofSample(models[[i]],re.form = NULL)  %>% 
    posterior_summary(probs=c(0.025,0.975),robust=T) %>% 
    as_tibble() 
  R2fix  <- R2_over_VarTot_OutofSample(models[[i]],re.form=NA)  %>% 
    posterior_summary(probs=c(0.025,0.975),robust=T) %>% 
    as_tibble() 
  df[names(models[i]),"R2all_Test"] <- paste0(round(R2all$Estimate, 3)," [",
                                              round(R2all$Q2.5, 3),"-",
                                              round(R2all$Q97.5, 3),"]") #  R2 of random + fixed effects
  df[names(models[i]),"R2fix_Test"] <- paste0(round(R2fix$Estimate, 3)," [",
                                              round(R2fix$Q2.5, 3),"-",
                                              round(R2fix$Q97.5, 3),"]")  # R2 of fixed effects
}


# Residuals and predictive errors    #####
# ======================================="
# Comments
# predictive_error uses the method "posterior_predict", so predict with considering sigma
# residuals uses the method "pp_expect", so give predicted means (without considering sigma)


for(m in 1:length(models)){
  
  # TRAIN - Predictive error
  ##########################"
  err_all <- predictive_error(models[[m]])  %>% 
    t() %>%  
    rowMeans() %>% 
    abs()  %>% 
    posterior_summary(probs=c(0.025,0.975)) %>% 
    as_tibble() 
  
  # Mean and CIs
  df[names(models[m]),"PE_Train"] <- paste0(round(err_all$Estimate,3), " [", 
                                            round(err_all$Q2.5,3),"-",
                                            round(err_all$Q97.5,3),"]")
  
  # # Mean and the standard error of the mean
  # df[names(models[m]),"PE_Train"] <- paste0(round(mean(err_all),3), " [", 
  #                                           round((sd(err_all)/sqrt(length(err_all))),3),"]")
  # 
  # # Mean and the standard deviation
  # df[names(models[m]),"PEsd_Train"] <- paste0(round(mean(err_all),3), " [", 
  #                                             round((sd(err_all)),3),"]")
  
  # TRAIN - Residuals
  ###################"
  # err_all <- residuals(models[[m]],method="pp_expect",summary=F)  %>% t() %>%  rowMeans() %>% abs()
  # 
  # # Mean and the standard error of the mean
  # df[names(models[m]),"RES_Train"] <- paste0(round(mean(err_all),3), " [", 
  #                                            round((sd(err_all)/sqrt(length(err_all))),3),"]")
  
  # TEST - Predictive error
  #########################"
  err_all <- predictive_error(models[[m]],newdata = test,allow_new_levels=T)  %>% 
    t() %>% 
    rowMeans() %>% 
    abs() %>% 
    posterior_summary(probs=c(0.025,0.975)) %>% 
    as_tibble() 
  
  # Mean and CIs
  df[names(models[m]),"PE_Test"] <- paste0(round(err_all$Estimate,3), " [", 
                                           round(err_all$Q2.5,3),"-",
                                           round(err_all$Q97.5,3),"]")
  
  # # Mean and the standard error of the mean
  # df[names(models[m]),"PE_Test"] <- paste0(round(mean(err_all),3), " [", 
  #                                          round((sd(err_all)/sqrt(length(err_all))),3),"]")
  # 
  # # Mean and the standard deviation
  # df[names(models[m]),"PEsd_Test"] <- paste0(round(mean(err_all),3), " [", 
  #                                            round((sd(err_all)),3),"]")
  # 
  # TEST - Residuals
  ##################"
  # err_all <- residuals(models[[m]],newdata = test,allow_new_levels=T,method="pp_expect",summary=F)  %>% t() %>%  rowMeans() %>% abs()
  # 
  # # Mean and the standard error of the mean
  # df[names(models[m]),"RES_Test"] <- paste0(round(mean(err_all),3), " [", 
  #                                           round((sd(err_all)/sqrt(length(err_all))),3),"]")
  
}


saveRDS(df, file=paste0("outputs/PerfTables/",part,"_ModelPerf_AfterReview.rds"))


# >>> latex tables for the manuscript ####
df <- readRDS(file=paste0("outputs/PerfTables/",part,"_ModelPerf_AfterReview.rds"))
df$Models <- paste0("M",str_sub(row.names(df),4,-1))

df <- df %>% 
  dplyr::select(Models,
                R2condAge_Train,
                R2all_Train,
                R2fix_Train,
                PE_Train,
                R2condAge_Test,
                R2all_Test,
                R2fix_Test,
                PE_Test)

# Generate the latex table
print(xtable(df, type = "latex",digits=0), file = paste0("tables/ModelPerf/",part,"_PerfModel_AfterReview.tex"), 
      include.rownames=FALSE)




# SITE-SPECIFIC TABLES (all sites merged) ####
# ============================================"
df <- readRDS(file=paste0("outputs/PerfTables/",part,"_ModelPerf_AfterReview.rds"))

# Extract the variance associated with age in each common garden:
var_age_baseline <- lapply(unique(train$site), function(site)
  pp_expect(models[[baseline]], transform = TRUE,newdata=test[test$site==site,],allow_new_levels=T,re.form=NA))
var_age_baseline <- lapply(var_age_baseline,function(x) apply(x,1,var)) 
var_age_baseline <-lapply(var_age_baseline,mean)
names(var_age_baseline) <- unique(train$site)
# Comment: the variance explained by age is null in Caceres and Madrid as there is only one age!


### Function to extract R2 without age in each site
R2_sites <- function(fit,site) {
  y <- log(test$height[test$site==site])
  var_y <- var(y)
  ypred <- pp_expect(fit, transform = TRUE,newdata=test[test$site==site,],allow_new_levels=T)
  var_ypred <- apply(ypred, 1, var)
  R2 <- (var_ypred - var_age_baseline[[site]]) / (var_y - var_age_baseline[[site]])
}


for (i in names(models)[!names(models) %in% c("MOD13")]){
  for(s in unique(models[[i]]$data$site)){
    
    R2  <- R2_sites(fit=models[[i]],site=s) %>% 
      posterior_summary(probs=c(0.025,0.975)) %>% 
      as_tibble() 
      
    df[names(models[i]),paste0("R2_Test_",s)] <- paste0(round(R2$Estimate,3), " [", 
                                                        round(R2$Q2.5,3),"-",
                                                        round(R2$Q97.5,3),"]")
  }}


df$Models <- paste0("M",str_sub(row.names(df),4,-1))

df <- df %>% 
  dplyr::select(Models,
                R2_Test_asturias,
                R2_Test_bordeaux,
                R2_Test_caceres,
                R2_Test_madrid,
                R2_Test_portugal)

# Generate the latex table
print(xtable(df, type = "latex",digits=0), file = paste0("tables/ModelPerf/",part,"_PerfModel_SiteSpecific_AfterReview.tex"), 
      include.rownames=FALSE)

