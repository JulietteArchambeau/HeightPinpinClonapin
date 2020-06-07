###############################################################################################"
##################                                                      #######################"
##################           Model performance - Figures                #######################"
##################                                                      #######################"
###############################################################################################"


library(loo)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(stringr)
library(xtable)
library(viridis)
library(ggplot2)
library(latex2exp)
library(brms)
library(gdata)

# my palette
if (part=="P1"|part=="P2")  mypalette <- c("#D73027", "#F46D43", "#FDAE61" ,"#A6DBA0" ,"#5AAE61","#E0F3F8" , "#ABD9E9", "#74ADD1", "#4575B4")
if (part=="P3") mypalette <- c("#D73027", "#F46D43", "#FDAE61" ,"#A6DBA0" ,"#5AAE61", "#74ADD1", "#4575B4")
  
# Data partition used to fit the models:
part <- "P3"

if (part=="P1"|part=="P2") list.models <- c(paste0("M",c(0:2,7:12)))
if (part=="P3") list.models <- c(paste0("M",c(0:2,7:8,11:12)))

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




# Site-specific figs ####

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

#saveRDS(df, file=paste0("outputs/PerfTables/",part,"_ModelPerf_SiteSpecific_R2_PE_GGplotTable.rds"))
df <- readRDS(file=paste0("outputs/PerfTables/",part,"_ModelPerf_SiteSpecific_R2_PE_GGplotTable.rds"))
df <- df[df$Models %in% list.models, ]

pR2 <- df %>% 
    mutate(Models = fct_relevel(Models, 
                                list.models)) %>% 
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
                              list.models)) %>% 
  ggplot(aes(fill=Models, y=PEmean, x=Sites)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  theme_bw() +
  geom_linerange(aes(ymin=PEmean-PEsd, ymax=PEmean+PEsd), position=position_dodge(.9),colour="gray30", alpha=0.7, size=1) + 
  labs(y=TeX("Site-specific mean predictive error"),x="") +
  scale_fill_manual(values=mypalette) +
  theme(legend.position = c(0.5,0.9),axis.text = element_text(size=22),axis.title = element_text(size=22), 
        legend.text = element_text(size=18),legend.title = element_text(size=20)) +
  guides(fill=guide_legend(ncol=3))


pPEse <- df %>% 
  mutate(Models = fct_relevel(Models, 
                              list.models)) %>% 
  ggplot(aes(fill=Models, y=PEmean, x=Sites)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  theme_bw() +
  geom_linerange(aes(ymin=PEmean-PEse, ymax=PEmean+PEse), position=position_dodge(.9),colour="gray30", alpha=0.7, size=1) + 
  labs(y=TeX("Site-specific mean predictive error"),x="") +
  scale_fill_manual(values=mypalette) +
  theme(legend.position = c(0.2,0.9),axis.text = element_text(size=22),axis.title = element_text(size=22), 
        legend.text = element_text(size=18),legend.title = element_text(size=20)) +
  guides(fill=guide_legend(ncol=3))


ggsave(pPEse,file=paste0("figs/SuppInfo/PEsecomparisons",part,".png"),height=8,width=15)


# > Figure in the manuscript ####

# with P1 and P2 partition

# P1 partition
df <- readRDS(file="outputs/PerfTables/P1_ModelPerf_SiteSpecific_R2_PE_GGplotTable.rds")
df <- df[df$Models %in% c(paste0("M",c(0:2,7:12))), ]

pR2P1 <- df %>% 
  mutate(Models = fct_relevel(Models, 
                              paste0("M",c(0:2,7:12)))) %>% 
  ggplot(aes(fill=Models, y=R2mean, x=Sites)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  theme_bw() +
  scale_fill_manual(values=mypalette) +
  geom_linerange(aes(ymin=R2mean-R2sd, ymax=R2mean+R2sd), position=position_dodge(.9),colour="gray30", alpha=0.7, size=1) + 
  labs(y=TeX("Site-specific $R^{2}_{s}$"),x="") +
  theme(legend.position = c(0.5,0.8),axis.text = element_text(size=22),axis.title = element_text(size=22), 
        legend.text = element_text(size=18),legend.title = element_text(size=20)) +
  guides(fill=guide_legend(ncol=3))

# P2 partition
df <- readRDS(file="outputs/PerfTables/P2_ModelPerf_SiteSpecific_R2_PE_GGplotTable.rds")

pR2P2 <- df %>% 
  mutate(Models = fct_relevel(Models, 
                              paste0("M",c(0:2,7:12)))) %>% 
  ggplot(aes(fill=Models, y=R2mean, x=Sites)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  theme_bw() +
  scale_fill_manual(values=mypalette) +
  geom_linerange(aes(ymin=R2mean-R2sd, ymax=R2mean+R2sd), position=position_dodge(.9),colour="gray30", alpha=0.7, size=1) + 
  labs(y=TeX("Site-specific $R^{2}_{s}$"),x="") +
  theme(legend.position = c("none"),axis.text = element_text(size=22),axis.title = element_text(size=22), 
        legend.text = element_text(size=18),legend.title = element_text(size=20)) +
  guides(fill=guide_legend(ncol=3))


fig <- ggarrange(pR2P1,pR2P2,labels=c("A)","B)"),font.label = list(size=20),nrow = 2)
ggsave(fig,file=paste0("figs/manuscript/R2comparisons.png"),height=10,width=15)



# Provenance-specific figs ####

if (grepl("P2",path)==TRUE|grepl("P3",path)==TRUE){
provs <- rep(unique(test$prov),length(names(models)))
vec <- c()
for(m in names(models)) vec <- c(vec,rep(m,6))
df <- data.frame(Models=vec,Provenances=provs,R2mean=NA,R2sd=NA)


for (i in 1:length(models)){
  
  for(p in unique(test$prov)){
    
    R2  <- brms::bayes_R2(models[[i]],re.form=NULL,newdata=test[test$prov==p,],allow_new_levels=T)
    df[df$Models==names(models[i])&df$Provenances==p,"R2mean"] <- R2[[1]]
    df[df$Models==names(models[i])&df$Provenances==p,"R2sd"] <- R2[[2]]
  }
}


for(i in 1:length(models)){
  
  test$pe <- predictive_error(models[[i]],newdata = test,allow_new_levels=T)  %>% t() %>%   rowMeans() %>% abs()
  
  for(p in unique(test$prov)){
    pe_prov <- test$pe[test$prov==p]
    
    df[df$Models==names(models[i])&df$Provenances==p,"PEmean"] <- mean(pe_prov)
    df[df$Models==names(models[i])&df$Provenances==p,"PEsd"] <- sd(pe_prov)
    df[df$Models==names(models[i])&df$Provenances==p,"PEse"] <- sd(pe_prov)/sqrt(length(pe_prov))
    
    
  }
  
}

df
df$Models <- paste0("M",str_sub(df$Models,4,-1))

saveRDS(df, file=paste0("outputs/PerfTables/",part,"_ModelPerf_ProvSpecific_R2_PE_GGplotTable.rds"))
df <- readRDS(file=paste0("outputs/PerfTables/",part,"_ModelPerf_ProvSpecific_R2_PE_GGplotTable.rds"))

df <- df[df$Models %in% c(paste0("M",c(0:2,7:12))), ]

pR2 <- df %>% 
  mutate(Models = fct_relevel(Models, 
                              list.models)) %>% 
  ggplot(aes(fill=Models, y=R2mean, x=Provenances)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  theme_bw() +
  scale_fill_manual(values=mypalette) +
  geom_linerange(aes(ymin=R2mean-R2sd, ymax=R2mean+R2sd), position=position_dodge(.9),colour="gray30", alpha=0.7, size=1) + 
  labs(y=TeX("Provenance-specific $R^{2}_{p}$"),x="") +
  theme(axis.text = element_text(size=22),axis.title = element_text(size=22), 
        legend.text = element_text(size=18),legend.title = element_text(size=20)) +
  guides(fill=guide_legend(ncol=1))

ggsave(pR2,file=paste0("figs/SuppInfo/R2comparisonsProv",part,".png"),height=8,width=15)

pPEsd <- df %>% 
  mutate(Models = fct_relevel(Models, 
                              list.models)) %>% 
  ggplot(aes(fill=Models, y=PEmean, x=Provenances)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  theme_bw() +
  geom_linerange(aes(ymin=PEmean-PEsd, ymax=PEmean+PEsd), position=position_dodge(.9),colour="gray30", alpha=0.7, size=1) + 
  labs(y=TeX("Provenance-specific $R^{2}$"),x="") +
  scale_fill_manual(values=mypalette) +
  theme(axis.text = element_text(size=22),axis.title = element_text(size=22), 
        legend.text = element_text(size=18),legend.title = element_text(size=20)) +
  guides(fill=guide_legend(ncol=1))


pPEse <- df %>% 
  mutate(Models = fct_relevel(Models, 
                              list.models)) %>% 
  ggplot(aes(fill=Models, y=PEmean, x=Provenances)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  theme_bw() +
  geom_linerange(aes(ymin=PEmean-PEse, ymax=PEmean+PEse), position=position_dodge(.9),colour="gray30", alpha=0.7, size=1) + 
  labs(y=TeX("Provenance-specific mean predictive error"),x="") +
  scale_fill_manual(values=mypalette) +
  theme(axis.text = element_text(size=22),axis.title = element_text(size=22), 
        legend.text = element_text(size=18),legend.title = element_text(size=20)) +
  guides(fill=guide_legend(ncol=1))


ggsave(pPEse,file=paste0("figs/SuppInfo/PEsecomparisonsProv",part,".png"),height=8,width=15)

}
