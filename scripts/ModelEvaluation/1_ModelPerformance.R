###############################################################################################"
##################                                                      #######################"
##################                  Model performance                   #######################"
##################                                                      #######################"
###############################################################################################"


library(loo)
library(dplyr)
library(tidyverse)
library(stringr)
library(xtable)
library(brms)
library(gdata)

# Data partition used to fit the models:
part <- "P2"

# Path
path= paste0("outputs/models/",part,"/")


# Models in a list

myFiles <- list.files(path=path,pattern=".rds")
models <- list()
for (i in 1:length(myFiles)){
  models[[i]] <- readRDS(file=paste0(path,myFiles[i]))
}
names(models) <- str_sub(myFiles,0,-5)



# test dataset

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



# 1) LOOIC ####

  # a) Compute looic only for new models ####
myloos <- list.files(paste0("outputs/loos/",part,"/"),pattern=".rds")
newmodels <- setdiff(str_sub(myFiles,1,-5),str_sub(myloos,5,-5)) 
new_loos <- list()
for (i in 1:length(newmodels)){
  new_loos[[i]] <- readRDS(file=paste0(path,newmodels[i],".rds"))
}

for (i in 1:length(new_loos)){
  myloo <- loo(new_loos[[i]],save_psis = F)
  saveRDS(myloo, file=paste0("outputs/loos/",part,"/loo_",newmodels[i],".rds"))
}


  # b) Looic table ####
myloos <- list.files(path=paste0("outputs/loos/",part,"/"),pattern=".rds") 
myloos <- str_sub(myloos,5,-5)
loos <- list()

for (i in 1:length(myloos)){
  loos[[i]] <- readRDS(file=paste0("outputs/loos/",part,"/","loo_",myloos[i],".rds"))
}
names(loos) <- myloos
df <- as.data.frame(loo::loo_compare(loos))
df



  # c) Looic pairwise comparisons ####

myloos <- list.files(path=paste0("outputs/loos/",part,"/"),pattern=".rds")
myloos <- str_sub(myloos,5,-5)
tab <- as.data.frame(matrix(NA,(length(myloos)-1),
                            (length(myloos)-1),
                            dimnames = list(c(myloos[1:length(myloos)-1]),c(myloos[2:length(myloos)]))))
loos <- list()
for (i in 1:length(myloos)){
  loos[[i]] <- readRDS(file=paste0("outputs/loos/",part,"/loo_",myloos[i],".rds"))
}
names(loos) <- myloos


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



if (grepl("P1",path)==TRUE){
colnames(tab) <- paste0("M",2:length(myloos))
rownames(tab) <- paste0("M",1:(length(myloos)-1))

} else if (grepl("P2",path)==TRUE){

colnames(tab) <- c("M7","M8","M2")
rownames(tab) <- c("M1","M7","M8")
}

print(xtable(tab, type = "latex",digits=2), file = paste0("tables/loos/",part,"_loo_pairwise_compare.tex"))




# 2) R2 ####

# Variance explained on the same data (train data set)
# 2 values are given: the mean of which is used as the measure of central tendency and the standard deviation as the measure of variability

for (i in 1:length(models)){
  R2all  <- brms::bayes_R2(models[[i]],re.form=NULL)
  R2fix  <- brms::bayes_R2(models[[i]],re.form=NA)
  df[names(models[i]),"R2all_Train"] <- paste0(round(R2all[[1]], 3)," [",round(R2all[[2]], 3),"]") #  R2 of random + fixed effects
  df[names(models[i]),"R2fix_Train"] <- paste0(round(R2fix[[1]], 3)," [",round(R2all[[2]], 3),"]")  # R2 of fixed effects
}


# Variance explained on new observations or new provenances (test data set)

for (i in 1:length(models)){
  R2all  <- brms::bayes_R2(models[[i]],re.form=NULL,newdata=test,allow_new_levels=T)
  R2fix  <- brms::bayes_R2(models[[i]],re.form=NA,newdata=test,allow_new_levels=T)
  df[names(models[i]),"R2all_Test"] <- paste0(round(R2all[[1]], 3)," [",
                                              round(R2all[[2]], 3),"]") # R2 of random + fixed effects
  df[names(models[i]),"R2fix_Test"] <- paste0(round(R2fix[[1]], 3)," [",
                                              round(R2all[[2]], 3),"]")  # R2 of fixed effects
}





# 3) Residuals and predictive errors #####

  # a) In all sites ####

# Comments
# predictive_error uses the method "posterior_predict", so predict with considering sigma
# residuals uses the method "pp_expect", so give predicted means (without considering sigma)


for(m in 1:length(models)){
  
  # TRAIN - Predictive error
  ##########################"
  err_all <- predictive_error(models[[m]])  %>% t() %>%  rowMeans() %>% abs()
  
  # Mean and the standard error of the mean
  df[names(models[m]),"PE_Train"] <- paste0(round(mean(err_all),3), " [", 
                                             round((sd(err_all)/sqrt(length(err_all))),3),"]")
  
  # TRAIN - Residuals
  ###################"
  err_all <- residuals(models[[m]],method="pp_expect",summary=F)  %>% t() %>%  rowMeans() %>% abs()
  
  # Mean and the standard error of the mean
  df[names(models[m]),"RES_Train"] <- paste0(round(mean(err_all),3), " [", 
                                             round((sd(err_all)/sqrt(length(err_all))),3),"]")
  
  # TEST - Predictive error
  #########################"
  err_all <- predictive_error(models[[m]],newdata = test,allow_new_levels=T)  %>% t() %>% rowMeans() %>% abs()
  
  # Mean and the standard error of the mean
  df[names(models[m]),"PE_Test"] <- paste0(round(mean(err_all),3), " [", 
                                             round((sd(err_all)/sqrt(length(err_all))),3),"]")
  
  # TEST - Residuals
  ##################"
  err_all <- residuals(models[[m]],newdata = test,allow_new_levels=T,method="pp_expect",summary=F)  %>% t() %>%  rowMeans() %>% abs()
  
  # Mean and the standard error of the mean
  df[names(models[m]),"RES_Test"] <- paste0(round(mean(err_all),3), " [", 
                                            round((sd(err_all)/sqrt(length(err_all))),3),"]")
  
}


saveRDS(df, file=paste0("outputs/PerfTables/",part,"_ModelPerf.rds"))





# >>> latex tables for the manuscript ####
df <- readRDS(file=paste0("outputs/PerfTables/",part,"_ModelPerf.rds"))
df$Models <- paste0("M",str_sub(row.names(df),4,-1))

df <- df %>% mutate(ELPD=paste0(round(elpd_loo,0)," [",round(se_elpd_loo,0),"]")) %>% 
    dplyr::select(Models,
                  ELPD=ELPD,
                  R2all_Train,
                  R2fix_Train,
                  PE_Train,
                  R2all_Test,
                  PE_Test)

# Generate the latex table
print(xtable(df, type = "latex",digits=0), file = paste0("tables/ModelPerf/",part,"_PerfModel.tex"), include.rownames=FALSE)



  # b) In each site ####
df <- readRDS(file=paste0("outputs/PerfTables/",part,"_ModelPerf.rds"))


for(m in 1:length(models)){
  
  # TRAIN - Predictive error
  ##########################"
  train$pe <- predictive_error(models[[m]])  %>% t() %>%  rowMeans() %>% abs()
  
  
  # TRAIN - Residuals
  ###################"
  train$res <- residuals(models[[m]],method="pp_expect",summary=F)  %>% t() %>%   rowMeans() %>% abs()
  
  
  # TEST - Predictive error
  #########################"
  test$pe <- predictive_error(models[[m]],newdata = test,allow_new_levels=T)  %>% t() %>%   rowMeans() %>% abs()
  
  # TEST - Residuals
  ##################"
  test$res <- residuals(models[[m]],newdata = test,allow_new_levels=T,method="pp_expect",summary=F)  %>% 
    t() %>%   rowMeans() %>% abs()
  
  
  for(i in unique(models[[m]]$data$site)){
    pe_site <- train$pe[train$site==i]
    pe_test_site <- test$pe[test$site==i]
    
    res_site <- train$res[train$site==i]
    res_test_site <- test$res[test$site==i]
    
    # Mean and the standard error of the mean
    df[names(models[m]),paste0("PE_Train_",i)] <- paste0(round(mean(pe_site),3), " [", 
                                                          round((sd(pe_site)/sqrt(length(pe_site))),3),"]")
    
    df[names(models[m]),paste0("PE_Test_",i)] <- paste0(round(mean(pe_test_site),3), " [", 
                                                         round((sd(pe_test_site)/sqrt(length(pe_test_site))),3),"]")

    df[names(models[m]),paste0("RES_Train_",i)] <- paste0(round(mean(res_site),3), " [", 
                                                          round((sd(res_site)/sqrt(length(res_site))),3),"]")
    
    df[names(models[m]),paste0("RES_Test_",i)] <- paste0(round(mean(res_test_site),3), " [", 
                                                         round((sd(res_test_site)/sqrt(length(res_test_site))),3),"]")
    
  }
  
}
df <- df[,1:36]
saveRDS(df, file=paste0("outputs/PerfTables/",part,"_ModelPerf_SiteSpecific.rds"))
#rownames(df)[rownames(df) == "MOD12"] <- "MOD8"
#rownames(df)[rownames(df) == "MOD11"] <- "MOD7"
df

# c) In each provenance (test dataset of P2 and P3) ####
if (grepl("P2",path)==TRUE|grepl("P3",path)==TRUE){
for(m in 1:length(models)){
  
  # TEST - Predictive error
  #########################"
  test$pe <- predictive_error(models[[m]],newdata = test,allow_new_levels=T)  %>% t() %>%  rowMeans() %>%  abs() 
  
  # TEST - Residuals
  ##################"
  test$res <- residuals(models[[m]],newdata = test,allow_new_levels=T,method="pp_expect",summary=F)  %>% 
    t() %>% rowMeans() %>%  abs() 
  
  for(i in unique(test$prov)){
    pe_test_prov <- test$pe[test$prov==i]
    res_test_prov <- test$res[test$prov==i]
    
    # Mean and the standard error of the mean
    df[names(models[m]),paste0("PE_Test_",i)] <- paste0(round(mean(pe_test_prov),3), " [", 
                                                         round((sd(pe_test_prov)/sqrt(length(pe_test_prov))),3),"]")
    
    df[names(models[m]),paste0("RES_Test_",i)] <- paste0(round(mean(res_test_prov),3), " [", 
                                                         round((sd(res_test_prov)/sqrt(length(res_test_prov))),3),"]")
    }
}
  
  saveRDS(df, file=paste0("outputs/PerfTables/",part,"_ModelPerf_SiteSpecificProvSpecific.rds"))
}



# >>> latex tables for the manuscript ####


# Site specific
df <- readRDS(file=paste0("outputs/PerfTables/",part,"_ModelPerf_SiteSpecific.rds"))
df$Models <- paste0("M",str_sub(row.names(df),4,-1))
df <- df %>% mutate(ELPD=paste0(round(elpd_loo,0)," [",round(se_elpd_loo,0),"]")) %>% 
  dplyr::select(Models,
                PE_Test_portugal,
                PE_Test_madrid,
                PE_Test_caceres,
                PE_Test_bordeaux,
                PE_Test_asturias)
print(xtable(df, type = "latex",digits=0), file = paste0("tables/ModelPerf/",part,"_PerfModel_SiteSpecific.tex"), include.rownames=FALSE)




df <- df %>% mutate(ELPD=paste0(round(elpd_loo,0)," [",round(se_elpd_loo,0),"]")) %>% 
  dplyr::select(Models,
                ELPD=ELPD,
                R2all_Train,
                R2fix_Train,
                PE_Train,
                RES_Train,
                R2all_Test,
                PE_Test,
                RES_Test)

df <- df %>% mutate(ELPD=paste0(round(elpd_loo,0)," [",round(se_elpd_loo,0),"]")) %>% 
  dplyr::select(Models,
                PE_Test_asturias,PE_Test_bordeaux,PE_Test_caceres,PE_Test_madrid,PE_Test_portugal)

# P2 partition - Prov specific
df <- readRDS(file=paste0("outputs/PerfTables/",part,"_ModelPerf_SiteSpecificProvSpecific.rds"))
df$Models <- paste0("M",str_sub(row.names(df),4,-1))
df <- df %>% mutate(ELPD=paste0(round(elpd_loo,0)," [",round(se_elpd_loo,0),"]")) %>% 
  dplyr::select(Models,
                PE_Test_CAD,
                PE_Test_COC,
                PE_Test_MIM,
                PE_Test_PLE,
                PE_Test_QUA,
                PE_Test_SAC)
print(xtable(df, type = "latex",digits=0), file = paste0("tables/ModelPerf/",part,"_PerfModel_ProvSpecific.tex"), include.rownames=FALSE)


# P3 partition - Prov specific
df <- readRDS(file=paste0("outputs/PerfTables/",part,"_ModelPerf_SiteSpecificProvSpecific.rds"))
df$Models <- paste0("M",str_sub(row.names(df),4,-1))
df <- df %>% mutate(ELPD=paste0(round(elpd_loo,0)," [",round(se_elpd_loo,0),"]")) %>% 
  dplyr::select(Models,
                PE_Test_TAM,
                PE_Test_SEG,
                PE_Test_PIA,
                PE_Test_ORI,
                PE_Test_OLO,
                PE_Test_LAM)
print(xtable(df, type = "latex",digits=0), file = paste0("tables/ModelPerf/",part,"_PerfModel_ProvSpecifi.tex"), include.rownames=FALSE)
