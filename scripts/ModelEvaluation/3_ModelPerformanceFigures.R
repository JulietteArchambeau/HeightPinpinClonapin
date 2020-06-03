###############################################################################################"
##################                                                      #######################"
##################           Model performance - Figures                #######################"
##################                                                      #######################"
###############################################################################################"


library(loo)
library(dplyr)
library(tidyverse)
library(stringr)
library(xtable)
library(viridis)
library(ggplot2)
library(latex2exp)
library(brms)
library(gdata)

# my palette
mypalette <- c("#D73027", "#F46D43", "#FDAE61" ,"#A6DBA0" ,"#5AAE61","#E0F3F8" , "#ABD9E9", "#74ADD1", "#4575B4")

# Data partition used to fit the models:
part <- "P1"

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


sites <- rep(unique(train$site),length(names(models)))
vec <- c()
for(m in names(models)) vec <- c(vec,rep(m,5))
df <- data.frame(Models=vec,Sites=sites,R2mean=NA,R2sd=NA)


for (i in 1:length(models)){
  if(names(models[i])=="MOD13") next
  
  for(s in unique(models[[i]]$data$site)){
    
    R2  <- brms::bayes_R2(models[[i]],re.form=NULL,newdata=test[test$site==s,],allow_new_levels=T)
    df[df$Models==names(models[i])&df$Sites==s,"R2mean"] <- R2[[1]]
    df[df$Models==names(models[i])&df$Sites==s,"R2sd"] <- R2[[2]]
  }
}


for(i in 1:length(models)){
  
  if(names(models[i])=="MOD13") next
  
  # TEST - Predictive error
  #########################"
  test$pe <- predictive_error(models[[i]],newdata = test,allow_new_levels=T)  %>% t() %>%   rowMeans() %>% abs()
  
  
  for(s in unique(models[[i]]$data$site)){
    pe_site <- test$pe[test$site==s]
    
    df[df$Models==names(models[i])&df$Sites==s,"PEmean"] <- mean(pe_site)
    df[df$Models==names(models[i])&df$Sites==s,"PEsd"] <- sd(pe_site)
    df[df$Models==names(models[i])&df$Sites==s,"PEse"] <- sd(pe_site)/sqrt(length(pe_site))
    
    
  }
  
}

df
df$Models <- paste0("M",str_sub(df$Models,4,-1))
df$Sites <- str_to_title(df$Sites)

saveRDS(df, file=paste0("outputs/PerfTables/",part,"_ModelPerf_SiteSpecific_R2_PE_GGplotTable.rds"))

df <- df[df$Models %in% c(paste0("M",c(0:2,7:12))), ]

pR2 <- df %>% 
    mutate(Models = fct_relevel(Models, 
                                paste0("M",c(0:2,7:12)))) %>% 
    ggplot(aes(fill=Models, y=R2mean, x=Sites)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    theme_bw() +
  scale_fill_manual(values=mypalette) +
    geom_linerange(aes(ymin=R2mean-R2sd, ymax=R2mean+R2sd), position=position_dodge(.9),colour="gray30", alpha=0.7, size=1) + 
  labs(y=TeX("Site-specific $R^{2}$"),x="") +
  theme(legend.position = c(0.5,0.9),axis.text = element_text(size=22),axis.title = element_text(size=22), 
        legend.text = element_text(size=18),legend.title = element_text(size=20)) +
  guides(fill=guide_legend(ncol=3))
  
ggsave(pR2,file=paste0("figs/manuscript/R2comparisons",part,".png"),height=9,width=14)

pPEsd <- df %>% 
  mutate(Models = fct_relevel(Models, 
                              paste0("M",c(0:2,7:12)))) %>% 
  ggplot(aes(fill=Models, y=PEmean, x=Sites)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  theme_bw() +
  geom_linerange(aes(ymin=PEmean-PEsd, ymax=PEmean+PEsd), position=position_dodge(.9),colour="gray30", alpha=0.7, size=1) + 
  labs(y=TeX("Site-specific $R^{2}$"),x="") +
  scale_fill_manual(values=mypalette) +
  theme(legend.position = c(0.5,0.9),axis.text = element_text(size=22),axis.title = element_text(size=22), 
        legend.text = element_text(size=18),legend.title = element_text(size=20)) +
  guides(fill=guide_legend(ncol=3))


pPEse <- df %>% 
  mutate(Models = fct_relevel(Models, 
                              paste0("M",c(0:2,7:12)))) %>% 
  ggplot(aes(fill=Models, y=PEmean, x=Sites)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  theme_bw() +
  geom_linerange(aes(ymin=PEmean-PEse, ymax=PEmean+PEse), position=position_dodge(.9),colour="gray30", alpha=0.7, size=1) + 
  labs(y=TeX("Site-specific $R^{2}$"),x="") +
  scale_fill_manual(values=mypalette) +
  theme(legend.position = c(0.5,0.9),axis.text = element_text(size=22),axis.title = element_text(size=22), 
        legend.text = element_text(size=18),legend.title = element_text(size=20)) +
  guides(fill=guide_legend(ncol=3))


ggsave(pPEse,file=paste0("figs/SuppInfo/PEsecomparisons",part,".png"),height=9,width=14)

