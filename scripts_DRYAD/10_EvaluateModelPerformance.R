########################################################################################################################"
#                                                                                                                      #
#                   Evalutating model goodness-of-fit and predictive ability                                           #
#                                                                                                                      #
#                                   Juliette Archambeau                                                                #
#                                       18/03/2022                                                                     #
#                                                                                                                      #
########################################################################################################################"



library(loo)       # CRAN v2.2.0
library(dplyr)     # CRAN v1.0.0
library(tidyverse) # CRAN v1.3.0
library(stringr)   # CRAN v1.4.0
library(xtable)    # CRAN v1.8-4
library(brms)      # CRAN v2.11.1
library(gdata)     # CRAN v2.18.0


# Baseline model:
# >> In the paper, we use M2 as baseline model as it has the highest predictive ability among the models relying only on the common garden design.
baseline <- "MOD2" 



# Data partition used to fit the models:
part <- "P2" # choose between P1, P2 and P3

# Path to the models:
path= paste0("outputs/models/",part,"/")


# Models in a list:
myFiles <- list.files(path=path,pattern=".rds")
models <- list()
for (i in 1:length(myFiles)){
  models[[i]] <- readRDS(file=paste0(path,myFiles[i]))
}
names(models) <- str_sub(myFiles,0,-5)



# Traning and test datasets:
if (grepl("P1",path)==TRUE){
  test  <- read_csv("data_DRYAD/TestP1prepared.csv")
  train <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv") %>%  dplyr::filter(P1=="train")
} else if (grepl("P2",path)==TRUE){
  test  <- read_csv("data_DRYAD/TestP2prepared.csv")
  train <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv") %>%  dplyr::filter(P2=="train")
} else if (grepl("P3",path)==TRUE){
  test  <- read_csv("data_DRYAD/TestP3prepared.csv")
  train <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv") %>%  dplyr::filter(P3=="train")
}



# 1) Model performance across all sites  ####
#    ==================================  "

# A) LOO    ####
#    ===    "

# Compute ELPD_loo:
myloos <- list.files(paste0("outputs/loos/",part,"/"),pattern=".rds")
newmodels <- setdiff(str_sub(myFiles,1,-5),str_sub(myloos,5,-5)) # allow calculating elpd_loo only for new models
new_loos <- list()
for (i in 1:length(newmodels)){
  new_loos[[i]] <- readRDS(file=paste0(path,newmodels[i],".rds"))
}

for (i in 1:length(new_loos)){
  myloo <- loo(new_loos[[i]],save_psis = F)
  saveRDS(myloo, file=paste0("outputs/loos/",part,"/loo_",newmodels[i],".rds"))
}


# Load loo tables:
myloos <- list.files(path=paste0("outputs/loos/",part,"/"),pattern=".rds") 
myloos <- str_sub(myloos,5,-5)
loos <- list()

for (i in 1:length(myloos)){
  loos[[i]] <- readRDS(file=paste0("outputs/loos/",part,"/","loo_",myloos[i],".rds"))
}
names(loos) <- myloos
df <- as.data.frame(loo::loo_compare(loos))
df

# >> Table S6. ####
# In this table: the mean (elpd_loo) and standard deviation (se_elpd_loo, within brackets in the table) of the ELPD_loo from model M0 to M12
print(xtable(df[!(row.names(df)=="MOD13"),c("elpd_loo", "se_elpd_loo")], type = "latex",digits=2), 
      file = paste0("tables/loos/",part,"_elpd_loo.tex"))


# >> Tables S7 and S10. ####
tab <- as.data.frame(matrix(NA,(length(myloos)-1),
                            (length(myloos)-1),
                            dimnames = list(c(myloos[1:length(myloos)-1]),c(myloos[2:length(myloos)]))))

for(i in myloos[1:length(myloos)-1]){
  for(j in myloos[2:length(myloos)]){
    m <- as.data.frame(loo_compare(list(modi=loos[[i]],modj=loos[[j]])))
    if(abs(loo_compare(loos[[i]],loos[[j]])[[2]])>(4*loo_compare(loos[[i]],loos[[j]])[[4]])){
      if (m["modi","elpd_diff"]==0){
        tab[i,j] <- paste0(round(-loo_compare(loos[[i]],loos[[j]])[[2]],2), " [",round(loo_compare(loos[[i]],loos[[j]])[[4]],2),"] *")
      } else if (m["modj","elpd_diff"]==0){
        tab[i,j] <- paste0(round(loo_compare(loos[[i]],loos[[j]])[[2]],2), " [",round(loo_compare(loos[[i]],loos[[j]])[[4]],2),"] *")
      }
    } else {
      if (m["modi","elpd_diff"]==0){
        tab[i,j] <- paste0(round(-loo_compare(loos[[i]],loos[[j]])[[2]],2), " [",round(loo_compare(loos[[i]],loos[[j]])[[4]],2),"]")
      } else if (m["modj","elpd_diff"]==0){
        tab[i,j] <- paste0(round(loo_compare(loos[[i]],loos[[j]])[[2]],2), " [",round(loo_compare(loos[[i]],loos[[j]])[[4]],2),"]")
      }
    }
    
    
  }
}

tab <- as.matrix(tab)
lowerTriangle(tab) <- NA
tab <- as.data.frame(tab)

colnames(tab) <- paste0("M",str_sub(colnames(tab),4,-1))
rownames(tab) <- paste0("M",str_sub(rownames(tab),4,-1))

tab
tab <- tab[!(rownames(tab)=="M13"),!(colnames(tab)=="M13")] # remove M13 in the P1 partition

# We save the tables (S7 for the P1 partition and S10 for the P2 partition):
print(xtable(tab, type = "latex",digits=2), file = paste0("tables/loos/",part,"_loo_pairwise_compare.tex"))




# B) R2 - Proportion of variance explained in the test and train datasets  #####
#    ====================================================================  "


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
  
}

# Saving intermediate file:
saveRDS(df, file=paste0("outputs/PerfTables/",part,"_ModelPerf_AfterReview.rds"))
df <- readRDS(file=paste0("outputs/PerfTables/",part,"_ModelPerf_AfterReview.rds"))


# >> Tables S4, S9 and S12. ####
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


print(xtable(df, type = "latex",digits=0), 
      file = paste0("tables/ModelPerf/",part,"_PerfModel_AfterReview.tex"), 
      include.rownames=FALSE)




#  2) Site-specific model performance  ####
#  ==================================  "
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

# >> Tables S8, S11 and S13. ####
print(xtable(df, type = "latex",digits=0), 
      file = paste0("tables/ModelPerf/",part,"_PerfModel_SiteSpecific_AfterReview.tex"), 
      include.rownames=FALSE)

