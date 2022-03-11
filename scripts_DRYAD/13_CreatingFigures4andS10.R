########################################################################################################################"
#                                                                                                                      #
#               Figures of the proportion of variance explained and predicted conditional on age                       #
#                                    Figures 4 to S10                                                                  #
#                                                                                                                      #
#                                   Juliette Archambeau                                                                #
#                                       11/03/2022                                                                     #
#                                                                                                                      #
########################################################################################################################"


library(dplyr)     # CRAN v1.0.0
library(tidyverse) # CRAN v1.3.0
library(ggpubr)    # CRAN v0.2.1
library(stringr)   # CRAN v1.4.0
library(latex2exp) # CRAN v0.4.0
library(ggplot2)   # CRAN v3.3.1
library(brms)      # CRAN v2.11.1
library(ggExtra)   # CRAN v0.9


# Baseline model
# >> In the paper, we use M2 as baseline model as it has the highest predictive ability among the models relying only on the common garden design.
baseline <- "MOD2"


# I) Proportion of predicted variance conditional on age  ####
# ======================================================  "

# >> Function to calculate the proportion of predicted variance conditional on age  ####
# ================================================================================  "

CalculateR2pred <- function(part){
  path= paste0("outputs/models/",part,"/")
  
  # Load the models:
  myFiles <- list.files(path=path,pattern=".rds")
  models <- list()
  for (i in 1:length(myFiles)){
    models[[i]] <- readRDS(file=paste0(path,myFiles[i]))
  }
  names(models) <- str_sub(myFiles,0,-5)
  
  models <- models[c(baseline,"MOD7","MOD8","MOD9","MOD10","MOD11","MOD12")]
  
  
  # Load the train and test datasets:
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
  
  # Site-specific R2
  # ================"
  
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
    var_ypred <- mean(apply(ypred, 1, var))
    R2 <- (var_ypred - var_age_baseline[[site]]) / (var_y - var_age_baseline[[site]])
  }
  
  
  df <- list()
  
  for (i in unique(train$site)) {
    df[[i]] <- lapply(models, R2_sites,site=i) %>%
      bind_cols(site=i)}
  
  df <- df %>% 
    bind_rows() %>% 
    dplyr::select(site,everything()) %>% 
    mutate(partition=part)
}


#  >> Apply the function on the three partitions    ####
#  =============================================    " 
df <- lapply(c("P1","P2","P3"),CalculateR2pred)



#  >> Build Figure 4    ####
#  =================    " 
p <- df %>% 
  bind_rows() %>% 
  mutate(site=replace(site, site=="caceres", "Cáceres")) %>% 
  pivot_longer(cols=contains("MOD"),names_to = "model",values_to="var") %>% 
  mutate(model = paste0("M",str_sub(.$model,4,-1)),
         model = fct_relevel(model,c(paste0("M",str_sub(baseline,4,-1)),paste0("M",c(7:12)))),
         site=str_to_title(site)) %>% 
  
  ggplot(aes(site,model,fill=var))+
  geom_tile(color= "white",size=0.1,width =0.8) + 
  scale_fill_gradient(low = "gray94", high = "forestgreen",limits=c(0, 0.4), breaks=seq(0,0.4,by=0.05)) +
  facet_grid(rows = vars(partition),scales="fixed", space = "fixed") +
  theme_minimal()+
  labs(fill=TeX("prediction $R^{2}_{ms}|age$")) +
  removeGrid() + 
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size=25,angle = 0),
        axis.title = element_blank(), 
        legend.text = element_text(size=18),
        legend.title = element_text(size=18))  +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 20,title.vjust = 3))

ggsave(p,file=paste0("figs/manuscript/ModelComparisonP1P2P3AfterReview2_",baseline,".png"),width = 10,height = 12)


# I) Proportion of explained variance conditional on age  ####
# ======================================================  "

# >> Function to calculate the proportion of explained variance conditional on age  ####
# ================================================================================  "

CalculateR2exp <- function(part){
  path= paste0("outputs/models/",part,"/")
  
  # Load the models:
  myFiles <- list.files(path=path,pattern=".rds")
  models <- list()
  for (i in 1:length(myFiles)){
    models[[i]] <- readRDS(file=paste0(path,myFiles[i]))
  }
  names(models) <- str_sub(myFiles,0,-5)
  models <- models[c(baseline,"MOD7","MOD8","MOD9","MOD10","MOD11","MOD12")]
  
  # Load the train and test datasets:
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
  
  # Site-specific R2
  # ================"
  
  # Extract the variance associated with age in each common garden:
  var_age_baseline <- lapply(unique(train$site), function(site)
    pp_expect(models[[baseline]], transform = TRUE,newdata=models[[baseline]]$data[models[[baseline]]$data$site==site,],re.form=NA))
  var_age_baseline <- lapply(var_age_baseline,function(x) apply(x,1,var)) 
  var_age_baseline <-lapply(var_age_baseline,mean)
  names(var_age_baseline) <- unique(train$site)
  # Comment: the variance explained by age is null in Caceres and Madrid as there is only one age!
  
  
  ### Function to extract R2 without age in each site
  R2_sites <- function(fit,site) {
    y <- fit$data$`log(height)`[fit$data$site==site]
    var_y <- var(y)
    ypred <- pp_expect(fit, transform = TRUE,newdata=fit$data[fit$data$site==site,])
    var_ypred <- mean(apply(ypred, 1, var))
    R2 <- (var_ypred - var_age_baseline[[site]]) / (var_y - var_age_baseline[[site]])
  }
  
  
  df <- list()
  
  for (i in unique(train$site)) {
    df[[i]] <- lapply(models, R2_sites,site=i) %>%
      bind_cols(site=i)}
  
  
  df <- df %>% 
    bind_rows() %>% 
    dplyr::select(site,everything()) %>% 
    mutate(partition=part)
}


#  >> Apply the function on the three partitions    ####
#  =============================================    " 
df <- lapply(c("P1","P2","P3"),CalculateR2exp)


#  >> Build Figure S10    ####
#  ===================    " 
p <- df %>% 
  bind_rows() %>% 
  mutate(site=replace(site, site=="caceres", "Cáceres")) %>% 
  pivot_longer(cols=contains("MOD"),names_to = "model",values_to="var") %>% 
  mutate(model = paste0("M",str_sub(.$model,4,-1)),
         model = fct_relevel(model,c(paste0("M",str_sub(baseline,4,-1)),paste0("M",c(7:12)))),
         site=str_to_title(site)) %>% 
  
  ggplot(aes(site,model,fill=var))+
  geom_tile(color= "white",size=0.1,width =0.8) + 
  scale_fill_gradient(low = "gray94", high = "forestgreen",limits=c(0, 0.4), breaks=seq(0,0.4,by=0.05)) +
  facet_grid(rows = vars(partition),scales="fixed", space = "fixed") +
  theme_minimal()+
  labs(fill=TeX("$R^{2}_{ms}|age$")) +
  removeGrid() + 
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size=25,angle = 0),
        axis.title = element_blank(), 
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))  +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 20,title.vjust = 3))

ggsave(p,file=paste0("figs/manuscript/ModelComparison_ExplicativePart_AfterReview2_",baseline,".png"),width = 10,height = 12)
