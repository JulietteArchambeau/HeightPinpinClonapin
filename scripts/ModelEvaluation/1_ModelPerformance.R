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


path= "outputs/models/P1/" # or path= "outputs/models/P2/"


# Models in a list

myFiles <- list.files(path=path,pattern=".rds")
models <- list()
for (i in 1:length(myFiles)){
  models[[i]] <- readRDS(file=paste0(path,myFiles[i]))
}
names(models) <- str_sub(myFiles,0,-5)



# test dataset

if (grepl("P1",path)==TRUE){
  test <- readRDS(file="data/datatest_pheno-clim-soil_20012020_Sample2_25percent_prepared.RDS")
  train <- readRDS(file="data/datatrain_pheno-clim-soil_20012020_Sample2_75percent.RDS")
} else if (grepl("P2",path)==TRUE){
  test <- readRDS(file="data/datatest_pheno-clim-soil_20012020_Sample2_6provs_prepared.RDS")
  train <- readRDS(file="data/datatrain_pheno-clim-soil_20012020_Sample4_28provs.RDS")
}


# 1) LOOIC ####

  # a) Compute looic only for new models ####

myloos <- list.files("outputs/loos/",pattern=".RDS")
newmodels <- setdiff(str_sub(myFiles,1,-5),str_sub(myloos,5,-5)) 
new_loos <- list()
for (i in 1:length(newmodels)){
  new_loos[[i]] <- readRDS(file=paste0(path,newmodels[i],".rds"))
}

for (i in 1:length(new_loos)){
  myloo <- loo(new_loos[[i]],save_psis = F)
  saveRDS(myloo, file=paste0("outputs/loos/loo_",newmodels[i],".rds"))
}


  # b) Looic table ####

myloos <- list.files(path=paste0("outputs/loos/"),pattern=".rds")
myloos <- str_sub(myloos,5,-5)
loos <- list()
for (i in 1:length(myloos)){
  loos[[i]] <- readRDS(file=paste0("outputs/loos/","loo_",myloos[i],".rds"))
}
names(loos) <- myloos
df <- as.data.frame(loo::loo_compare(loos))
df



  # c) Looic pairwise comparisons ####

myloos <- list.files(path=paste0("outputs/loos/"),pattern=".rds")
myloos <- str_sub(myloos,5,-5)
tab <- as.data.frame(matrix(NA,(length(myloos)-1),
                            (length(myloos)-1),
                            dimnames = list(c(myloos[1:length(myloos)-1]),c(myloos[2:length(myloos)]))))
loos <- list()
for (i in 1:length(myloos)){
  loos[[i]] <- readRDS(file=paste0("outputs/loos/","loo_",myloos[i],".rds"))
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

print(xtable(tab, type = "latex",digits=2), file = "tables/loos/P1_loo_pairwise_compare.tex")

} else if (grepl("P2",path)==TRUE){

colnames(tab) <- c("M2","M3","M4","M6","M7","M8","M9","M10")
rownames(tab) <- c("M1","M2","M3","M4","M6","M7","M8","M9")

print(xtable(tab, type = "latex",digits=2), file = "tables/loos/P2_loo_pairwise_compare.tex")
}







# 2) R2 ####

# Variance explained on the same data (train data set)

for (i in 1:length(models)){
  R2all  <- brms::bayes_R2(models[[i]],re.form=NULL)
  R2fix  <- brms::bayes_R2(models[[i]],re.form=NA)
  df[names(models[i]),"R2all"] <- paste0(round(R2all[[1]], 3)," [",round(R2all[[2]], 3),"]") #  R2 of random + fixed effects
  df[names(models[i]),"R2fix"] <- paste0(round(R2fix[[1]], 3)," [",round(R2all[[2]], 3),"]")  # R2 of fixed effects
}


# Variance explained on new observations or new provenances (test data set)

for (i in 1:length(models)){
  R2all  <- brms::bayes_R2(models[[i]],re.form=NULL,newdata=test,allow_new_levels=T)
  R2fix  <- brms::bayes_R2(models[[i]],re.form=NA,newdata=test,allow_new_levels=T)
  df[names(models[i]),"R2all_test"] <- paste0(round(R2all[[1]], 3)," [",
                                              round(R2all[[2]], 3),"]") # R2 of random + fixed effects
  df[names(models[i]),"R2fix_test"] <- paste0(round(R2fix[[1]], 3)," [",
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
  err_all <- predictive_error(models[[m]])
  err_all <- tibble(err_all)  %>% t() %>%  abs(.) %>% rowMeans(.) %>% as.numeric()
  
  # Mean and the standard error of the mean
  df[names(models[m]),"PE_Train1"] <- paste0(round(mean(err_all),3), " [", 
                                             round((sd(err_all)/sqrt(length(err_all))),3),"]")
  
  # TRAIN - Residuals
  ###################"
  err_all <- residuals(models[[m]],method="pp_expect",summary=F)
  err_all <- tibble(err_all)  %>% t() %>%  abs(.) %>% rowMeans(.) %>% as.numeric()
  
  # Mean and the standard error of the mean
  df[names(models[m]),"RES_Train1"] <- paste0(round(mean(err_all),3), " [", 
                                             round((sd(err_all)/sqrt(length(err_all))),3),"]")
  
  # TEST - Predictive error
  #########################"
  err_all <- predictive_error(models[[m]],newdata = test,allow_new_levels=T)
  err_all <- tibble(err_all)  %>% t() %>%  abs(.) %>% rowMeans(.) %>% as.numeric()
  
  # Mean and the standard error of the mean
  df[names(models[m]),"PE_Test1"] <- paste0(round(mean(err_all),3), " [", 
                                             round((sd(err_all)/sqrt(length(err_all))),3),"]")
  
  # TEST - Residuals
  ##################"
  err_all <- residuals(models[[m]],newdata = test,allow_new_levels=T,method="pp_expect",summary=F)
  err_all <- tibble(err_all)  %>% t() %>%  abs(.) %>% rowMeans(.) %>% as.numeric()
  
  # Mean and the standard error of the mean
  df[names(models[m]),"RES_Test1"] <- paste0(round(mean(err_all),3), " [", 
                                            round((sd(err_all)/sqrt(length(err_all))),3),"]")
  
}



if (grepl("P1",path)==TRUE){
  
  saveRDS(df, file=paste0("outputs/PerfTables/P1_ModelPerf.rds"))
  
} else if (grepl("P2",path)==TRUE){
  
  saveRDS(df, file=paste0("outputs/PerfTables/P2_ModelPerf.rds"))
  
}




# >>> latex tables for the manuscript ####

df$Models <- paste0("M",1:length(models))

df <- df %>% mutate(ELPD=paste0(round(elpd_loo,0)," [",round(se_elpd_loo,0),"]")) %>% 
  dplyr::select(Models,
                ELPD=ELPD,
                R2all_train=R2all,
                R2fix_train=R2fix,
                PE_train=PE_Train1,
                R2all_test=R2all_test,
                PE_test=PE_Test1)


# Generate the latex table
if (grepl("P1",path)==TRUE){
  
  print(xtable(df, type = "latex",digits=0), file = "tables/ModelPerf/P1_PerfModel.tex", include.rownames=FALSE)
  
} else if (grepl("P2",path)==TRUE){
  
  print(xtable(df, type = "latex",digits=0), file = "tables/ModelPerf/P2_PerfModel.tex", include.rownames=FALSE)
  
}



  # b) In each site ####
if (grepl("P1",path)==TRUE){
  df <- readRDS(file="outputs/PerfTables/P1_ModelPerf.rds")
  
} else if (grepl("P2",path)==TRUE){
  
  df <- readRDS(file="outputs/PerfTables/P2_ModelPerf.rds")
}


for(m in 1:length(models)){
  
  # TRAIN - Predictive error
  ##########################"
  pe <- predictive_error(models[[m]])
  pe <- tibble(pe)  %>% t() %>%  abs(.) %>% rowMeans(.) %>% as.numeric()
  train$pe <- pe
  
  
  # TRAIN - Residuals
  ###################"
  res <- residuals(models[[m]],method="pp_expect",summary=F)
  res <- tibble(res)  %>% t() %>%  abs(.) %>% rowMeans(.) %>% as.numeric()
  train$res <- res
  
  
  # TEST - Predictive error
  #########################"
  pe_test <- predictive_error(models[[m]],newdata = test,allow_new_levels=T)
  pe_test <- tibble(pe_test)  %>% t() %>%  abs(.) %>% rowMeans(.) %>% as.numeric()
  test$pe <- pe_test
  
  # TEST - Residuals
  ##################"
  res_test <- residuals(models[[m]],newdata = test,allow_new_levels=T,method="pp_expect",summary=F)
  res_test <- tibble(res_test)  %>% t() %>%  abs(.) %>% rowMeans(.) %>% as.numeric()
  test$res <- res_test
  
  
  for(i in unique(models[[m]]$data$site)){
    pe_site <- train$pe[train$site==i]
    pe_test_site <- test$pe[test$site==i]
    
    res_site <- train$res[train$site==i]
    res_test_site <- test$res[test$site==i]
    
    # Mean and the standard error of the mean
    df[names(models[m]),paste0("PE_Train1_",i)] <- paste0(round(mean(pe_site),3), " [", 
                                                          round((sd(pe_site)/sqrt(length(pe_site))),3),"]")
    
    df[names(models[m]),paste0("PE_Test1_",i)] <- paste0(round(mean(pe_test_site),3), " [", 
                                                         round((sd(pe_test_site)/sqrt(length(pe_test_site))),3),"]")

    df[names(models[m]),paste0("RES_Train1_",i)] <- paste0(round(mean(res_site),3), " [", 
                                                          round((sd(res_site)/sqrt(length(res_site))),3),"]")
    
    df[names(models[m]),paste0("RES_Test1_",i)] <- paste0(round(mean(res_test_site),3), " [", 
                                                         round((sd(res_test_site)/sqrt(length(res_test_site))),3),"]")
    
  }
  
}



# c) In each provenance (test dataset of P2) ####
if (grepl("P2",path)==TRUE){
for(m in 1:length(models)){
  
  # TEST - Predictive error
  #########################"
  pe_test <- predictive_error(models[[m]],newdata = test,allow_new_levels=T)
  pe_test <- tibble(pe_test)  %>% t() %>%  abs(.) %>% rowMeans(.) %>% as.numeric()
  test$pe <- pe_test
  
  for(i in unique(test$prov)){
    pe_test_prov <- test$pe[test$prov==i]
    
    # Mean and the standard error of the mean
    df[names(models[m]),paste0("PE_Test1_",i)] <- paste0(round(mean(pe_test_prov),3), " [", 
                                                         round((sd(pe_test_prov)/sqrt(length(pe_test_prov))),3),"]")
    }
  }
}

if (grepl("P1",path)==TRUE){
  
  saveRDS(df, file=paste0("outputs/PerfTables/P1_ModelPerf_SiteSpecific.rds"))
  
} else if (grepl("P2",path)==TRUE){
  
  saveRDS(df, file=paste0("outputs/PerfTables/P2_ModelPerf_SiteProvSpecific.rds"))
  
}
