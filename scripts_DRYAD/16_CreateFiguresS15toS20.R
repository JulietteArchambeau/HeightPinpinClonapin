########################################################################################################################"
#                                                                                                                      #
#               Figures of the posterior distributions in models M7 to M12                                             #
#                                   Figures S15 to S20                                                                 #
#                                                                                                                      #
#                                   Juliette Archambeau                                                                #
#                                       18/03/2022                                                                     #
#                                                                                                                      #
########################################################################################################################"


library(brms)      # CRAN v2.11.1
library(ggplot2)   # CRAN v3.3.1
library(ggridges)  # CRAN v0.5.1
library(latex2exp) # CRAN v0.4.0
library(tibble)    # CRAN v2.1.3
library(ggpubr)    # CRAN v0.2.1
library(tidyverse) # CRAN v1.3.0



# Select the partition:
part <- "P1" # among P1, P2 and P3



# Function used in the script:
source("scripts/Functions/vir_lite.R") # available here: https://github.com/JulietteArchambeau/HeightPinpinClonapin/blob/master/scripts/Functions/vir_lite.R
ds=0.7 # parameter of the function vir_lite



# A) Posterior distributions in M7 and M8    #####
# =======================================    "


# Load models:
models <- list()
models[[1]] <- readRDS(file=paste0("outputs/models/",part,"/MOD7.rds"))
models[[2]] <- readRDS(file=paste0("outputs/models/",part,"/MOD8.rds"))
names(models) <- c("M7","M8")


# create lists to store figures
figs.GP <- list()
figs.beta <- list()


# >> Panels a: Population structure  ####
# ---------------------------------  "

for (i in 1:2){
  POST <- posterior_samples(models[[i]],pars = "^r_mmQ1Q2Q3Q4Q5Q6\\[")
  colnames(POST) <- str_sub(colnames(POST),18,-12)
  POST <- as.data.frame(t(POST))
  POST$genepool <- as.factor(rownames(POST))
  
  posteriorsimpelmodellong <- POST %>%  as_tibble() %>% 
    gather(key = "key", value = "value", -genepool)%>%
    group_by(genepool) %>%
    dplyr::mutate(meanpergenepool = mean(value))%>%
    ungroup()
  
  figs.GP[[i]] <- ggplot()+
    geom_vline(xintercept = 0, 
               col        = "grey70") +
    stat_density_ridges(data  = posteriorsimpelmodellong, 
                        aes(x      = value,
                            y      = reorder(as.factor(genepool), meanpergenepool),
                            fill   = as.factor(genepool),
                            vline_color = ..quantile..),
                        scale = 2, 
                        alpha = .6,
                        rel_min_height=c(.0044),
                        size=0.5,
                        quantile_lines = TRUE, quantiles = c(0.025,0.5,0.975)) +
    scale_discrete_manual("vline_color",
                          values = c("blue", "red", "blue", "black"), 
                          breaks = c(1, 2),
                          labels = c("2.5 and 97.5th percentiles", "Median"),
                          name = NULL) +
    scale_y_discrete(labels=c("Q1"=parse(text = TeX("$g_{NA}$")),
                              "Q2"=parse(text = TeX("$g_{C}$")),
                              "Q3"=parse(text = TeX("$g_{CS}$")),
                              "Q4"=parse(text = TeX("$g_{FA}$")),
                              "Q5"=parse(text = TeX("$g_{IA}$")),
                              "Q6"=parse(text = TeX("$g_{SES}$")))) +
    coord_cartesian(c(-0.5,0.6))+
    scale_fill_manual(values=c("orangered3",
                               "gold2",
                               "darkorchid3",
                               "navyblue",
                               "turquoise2",
                               "green3"), labels = c("Northern Africa (NA)", 
                                                     "Corsica (C)",
                                                     "Central Spain (CS)",
                                                     "French Atlantic (FA)",
                                                     "Iberian Atlantic (IA)",
                                                     "South-eastern Spain (SES)")) +
    labs(fill     = "Gene pools",
         y        = "", 
         x = ""
    ) + 
    theme_bw() + theme(axis.text = element_text(size=22),axis.title = element_text(size=22), 
                       legend.text = element_text(size=16),legend.title = element_text(size=17))
}

figs.GP[[1]] <- figs.GP[[1]] + theme(legend.position = "none")
pGP <- ggarrange(figs.GP[[1]],figs.GP[[2]],
                 labels=c("M7. a)","M8. a)"),
                 font.label = list(size = 20),
                 hjust=c(-0.1,-0.1),
                 vjust=c(1.6,1.6),
                 nrow=1,
                 widths = c(1,1.3))




# >> Panels b: Provenance climates and PEAs   ####
# -----------------------------------------   "


for(i in 1:2){
  
  if(i==1){
    pea="gPEA"
    variable ="gPEA.sc"
  } else if (i==2){
    pea="rPEA"
    variable ="rPEA.sc"
  }
  
  
  POST <- posterior_samples(models[[i]],pars = "^r_site\\[.*sc\\]") %>% dplyr::rename(
    beta_PEA_Portugal = paste0('r_site[portugal,',variable,']'),
    beta_PEA_Bordeaux = paste0('r_site[bordeaux,',variable,']'),
    beta_PEA_Asturias = paste0('r_site[asturias,',variable,']'),
    beta_PEA_Madrid = paste0('r_site[madrid,',variable,']'),
    beta_PEA_Caceres = paste0('r_site[caceres,',variable,']'),
    beta_MinPre_Portugal = 'r_site[portugal,bio14_prov.sc]',
    beta_MinPre_Bordeaux = 'r_site[bordeaux,bio14_prov.sc]',
    beta_MinPre_Asturias = 'r_site[asturias,bio14_prov.sc]',
    beta_MinPre_Madrid = 'r_site[madrid,bio14_prov.sc]',
    beta_MinPre_Caceres = 'r_site[caceres,bio14_prov.sc]',
    beta_MaxTemp_Portugal = 'r_site[portugal,bio5_prov.sc]',
    beta_MaxTemp_Bordeaux = 'r_site[bordeaux,bio5_prov.sc]',
    beta_MaxTemp_Asturias = 'r_site[asturias,bio5_prov.sc]',
    beta_MaxTemp_Madrid = 'r_site[madrid,bio5_prov.sc]',
    beta_MaxTemp_Caceres = 'r_site[caceres,bio5_prov.sc]'
  )
  POST <- as.data.frame(t(POST))
  POST$var <- as.factor(rownames(POST))
  
  posteriorsimpelmodellong <- POST %>%  as_tibble() %>% 
    gather(key = "key", value = "value", -var)
  
  
  figs.beta[[i]] <- ggplot()+
    geom_vline(xintercept = 0, col="grey70") +
    stat_density_ridges(data  = posteriorsimpelmodellong, 
                        aes(x      = value,
                            y      = var,
                            fill   = as.factor(var),
                            vline_color = ..quantile..),
                        scale = 1.8, 
                        alpha = .8,
                        size=0.5,
                        rel_min_height=.01,
                        quantile_lines = TRUE, 
                        quantiles = c(0.025,0.5,0.975)) +
    # coord_cartesian(c(-,0.8))+
    scale_discrete_manual("vline_color",
                          values = c("blue", "red", "blue", "black"), 
                          breaks = c(1, 2),
                          labels = c("5% & 95% quantiles", "mean"),
                          name = NULL) +
    scale_y_discrete(labels=c("beta_MaxTemp_Caceres"=parse(text = TeX("$\\beta_{max.temp,Caceres}$")),
                              "beta_MaxTemp_Bordeaux"=parse(text = TeX("$\\beta_{max.temp,Bordeaux}$")),
                              "beta_MaxTemp_Portugal"=parse(text = TeX("$\\beta_{max.temp,Portugal}$")),
                              "beta_MaxTemp_Madrid"=parse(text = TeX("$\\beta_{max.temp,Madrid}$")),
                              "beta_MaxTemp_Asturias"=parse(text = TeX("$\\beta_{max.temp,Asturias}$")),
                              "beta_MinPre_Caceres"=parse(text = TeX("$\\beta_{min.pre,Caceres}$")),
                              "beta_MinPre_Bordeaux"=parse(text = TeX("$\\beta_{min.pre,Bordeaux}$")),
                              "beta_MinPre_Portugal"=parse(text = TeX("$\\beta_{min.pre,Portugal}$")),
                              "beta_MinPre_Madrid"=parse(text = TeX("$\\beta_{min.pre,Madrid}$")),
                              "beta_MinPre_Asturias"=parse(text = TeX("$\\beta_{min.pre,Asturias}$")),
                              "beta_PEA_Caceres"=parse(text = TeX(paste0("$\\beta_{",pea,",Caceres}$"))),
                              "beta_PEA_Madrid"=parse(text = TeX(paste0("$\\beta_{",pea,",Madrid}$"))),
                              "beta_PEA_Portugal"=parse(text = TeX(paste0("$\\beta_{",pea,",Portugal}$"))),
                              "beta_PEA_Asturias"=parse(text = TeX(paste0("$\\beta_{",pea,",Asturias}$"))),
                              "beta_PEA_Bordeaux"=parse(text = TeX(paste0("$\\beta_{",pea,",Bordeaux}$")))
    )) + 
    labs(y = "", 
         x = ""
         #x        = TeX(paste0("Estimate of the varying slopes $\\beta_{",pea,",site}$, $\\beta_{min.pre,site}$ and $\\beta_{max.temp,site}$"))
    ) + 
    scale_fill_manual(values=c(vir_lite("cyan2",ds=ds),
                               vir_lite("navyblue",ds=ds),
                               vir_lite("pink",ds=ds),
                               vir_lite("deeppink",ds=ds),
                               vir_lite("dodgerblue2",ds=ds),
                               vir_lite("cyan2",ds=ds),
                               vir_lite("navyblue",ds=ds),
                               vir_lite("pink",ds=ds),
                               vir_lite("deeppink",ds=ds),
                               vir_lite("dodgerblue2",ds=ds),
                               vir_lite("cyan2",ds=ds),
                               vir_lite("navyblue",ds=ds),
                               vir_lite("pink",ds=ds),
                               vir_lite("deeppink",ds=ds),
                               vir_lite("dodgerblue2",ds=ds))) +
    theme_bw() + theme(axis.text = element_text(size=22),axis.title = element_text(size=22),
                       legend.position = "none",
                       legend.text = element_text(size=16),legend.title = element_text(size=17))
}

pbeta <- ggarrange(figs.beta[[1]],figs.beta[[2]],labels=c("M7. b)","M8. b)"),font.label = list(size = 20),hjust=c(-0.1,-0.1),vjust=c(1.6,1.6),nrow=1)


# >> Figures S15, S17 and S19 ####

figtot <- ggarrange(pGP,pbeta,nrow=2,heights=c(1,2)) 
figtot
ggsave(figtot,file=paste0("figs/SuppInfo/M7M8Posteriors",part,".png"),height=12,width=20)








#  B)  Posterior distributions in M9 to M12  ####
#  ========================================  "

# >>  M9 => Pop structure  ####
# -----------------------  "

mod <- readRDS(file=paste0("outputs/models/",part,"/MOD9.rds"))

POST <- posterior_samples(mod,pars = "^r_mmQ1Q2Q3Q4Q5Q6\\[")
colnames(POST) <- str_sub(colnames(POST),18,-12)
POST <- as.data.frame(t(POST))
POST$genepool <- as.factor(rownames(POST))

posteriorsimpelmodellong <- POST %>%  as_tibble() %>% 
  gather(key = "key", value = "value", -genepool)%>%
  group_by(genepool) %>%
  dplyr::mutate(meanpergenepool = mean(value))%>%
  ungroup()

pGP <- ggplot()+
  geom_vline(xintercept = 0, 
             col        = "grey70") +
  stat_density_ridges(data  = posteriorsimpelmodellong, 
                      aes(x      = value,
                          y      = reorder(as.factor(genepool), meanpergenepool),
                          fill   = as.factor(genepool),
                          vline_color = ..quantile..),
                      scale = 2, 
                      alpha = .6,
                      rel_min_height=c(.001),
                      size=0.5,
                      quantile_lines = TRUE, quantiles = c(0.025,0.5,0.975)) +
  scale_discrete_manual("vline_color",
                        values = c("blue", "red", "blue", "black"), 
                        breaks = c(1, 2),
                        labels = c("2.5 and 97.5th percentiles", "Median"),
                        name = NULL) +
  coord_cartesian(c(-0.5,0.4))+
  scale_fill_manual(values=c("orangered3",
                             "gold2",
                             "darkorchid3",
                             "navyblue",
                             "turquoise2",
                             "green3"), labels = c("Northern Africa (NA)", 
                                                   "Corsica (C)",
                                                   "Central Spain (CS)",
                                                   "French Atlantic (FA)",
                                                   "Iberian Atlantic (IA)",
                                                   "South-eastern Spain (SES)")) +
  scale_y_discrete(labels=c("Q1"=parse(text = TeX("$g_{NA}$")),
                            "Q2"=parse(text = TeX("$g_{C}$")),
                            "Q3"=parse(text = TeX("$g_{CS}$")),
                            "Q4"=parse(text = TeX("$g_{FA}$")),
                            "Q5"=parse(text = TeX("$g_{IA}$")),
                            "Q6"=parse(text = TeX("$g_{SES}$")))) + 
  labs(fill     = "Gene pools",
       y        = "", 
       x=""
  ) + 
  theme_bw() + theme(axis.text = element_text(size=18),axis.title = element_text(size=18), 
                     legend.text = element_text(size=16),legend.title = element_text(size=17))



# >> M10 => Climate in the provenances   ####
# -------------------------------------  "


mod <- readRDS(file=paste0("outputs/models/",part,"/MOD10.rds"))

POST <- posterior_samples(mod,pars = "^r_site\\[.*sc\\]") %>% dplyr::rename(
  beta_MinPre_Portugal = 'r_site[portugal,bio14_prov.sc]',
  beta_MinPre_Bordeaux = 'r_site[bordeaux,bio14_prov.sc]',
  beta_MinPre_Asturias = 'r_site[asturias,bio14_prov.sc]',
  beta_MinPre_Madrid = 'r_site[madrid,bio14_prov.sc]',
  beta_MinPre_Caceres = 'r_site[caceres,bio14_prov.sc]',
  beta_MaxTemp_Portugal = 'r_site[portugal,bio5_prov.sc]',
  beta_MaxTemp_Bordeaux = 'r_site[bordeaux,bio5_prov.sc]',
  beta_MaxTemp_Asturias = 'r_site[asturias,bio5_prov.sc]',
  beta_MaxTemp_Madrid = 'r_site[madrid,bio5_prov.sc]',
  beta_MaxTemp_Caceres = 'r_site[caceres,bio5_prov.sc]'
)
POST <- as.data.frame(t(POST))
POST$var <- as.factor(rownames(POST))

posteriorsimpelmodellong <- POST %>%  as_tibble() %>% 
  gather(key = "key", value = "value", -var)


pCP <- ggplot()+
  geom_vline(xintercept = 0, col="grey70") +
  stat_density_ridges(data  = posteriorsimpelmodellong, 
                      aes(x      = value,
                          y      = var,
                          fill   = as.factor(var),
                          vline_color = ..quantile..),
                      scale = 1.8, 
                      alpha = .8,
                      size=0.5,
                      rel_min_height=.01,
                      quantile_lines = TRUE, 
                      quantiles = c(0.025,0.5,0.975)) +
  scale_discrete_manual("vline_color",
                        values = c("blue", "red", "blue", "black"), 
                        breaks = c(1, 2),
                        labels = c("5% & 95% quantiles", "mean"),
                        name = NULL) +
  scale_y_discrete(labels=c("beta_MaxTemp_Caceres"=parse(text = TeX("$\\beta_{max.temp,Caceres}$")),
                            "beta_MaxTemp_Bordeaux"=parse(text = TeX("$\\beta_{max.temp,Bordeaux}$")),
                            "beta_MaxTemp_Portugal"=parse(text = TeX("$\\beta_{max.temp,Portugal}$")),
                            "beta_MaxTemp_Madrid"=parse(text = TeX("$\\beta_{max.temp,Madrid}$")),
                            "beta_MaxTemp_Asturias"=parse(text = TeX("$\\beta_{max.temp,Asturias}$")),
                            "beta_MinPre_Caceres"=parse(text = TeX("$\\beta_{min.pre,Caceres}$")),
                            "beta_MinPre_Bordeaux"=parse(text = TeX("$\\beta_{min.pre,Bordeaux}$")),
                            "beta_MinPre_Portugal"=parse(text = TeX("$\\beta_{min.pre,Portugal}$")),
                            "beta_MinPre_Madrid"=parse(text = TeX("$\\beta_{min.pre,Madrid}$")),
                            "beta_MinPre_Asturias"=parse(text = TeX("$\\beta_{min.pre,Asturias}$")))) + 
  labs(y        = "", 
       x="") +
  scale_fill_manual(values=c(vir_lite("cyan2",ds=ds),
                             vir_lite("navyblue",ds=ds),
                             vir_lite("pink",ds=ds),
                             vir_lite("deeppink",ds=ds),
                             vir_lite("dodgerblue2",ds=ds),
                             vir_lite("cyan2",ds=ds),
                             vir_lite("navyblue",ds=ds),
                             vir_lite("pink",ds=ds),
                             vir_lite("deeppink",ds=ds),
                             vir_lite("dodgerblue2",ds=ds))) +
  theme_bw() + theme(axis.text = element_text(size=22),axis.title = element_text(size=22), 
                     legend.position = "none", plot.title = element_text(size=22)) 



# >> M11 => gPEAs  ####
# ---------------- "

mod <- readRDS(file=paste0("outputs/models/",part,"/MOD11.rds"))

pea="gPEA"
variable ="gPEA.sc"

POST <- posterior_samples(mod,pars = "^r_site\\[.*sc\\]") %>% dplyr::rename(
  beta_PEA_Portugal = paste0('r_site[portugal,',variable,']'),
  beta_PEA_Bordeaux = paste0('r_site[bordeaux,',variable,']'),
  beta_PEA_Asturias = paste0('r_site[asturias,',variable,']'),
  beta_PEA_Madrid = paste0('r_site[madrid,',variable,']'),
  beta_PEA_Caceres = paste0('r_site[caceres,',variable,']')
)
POST <- as.data.frame(t(POST))
POST$var <- as.factor(rownames(POST))

posteriorsimpelmodellong <- POST %>%  as_tibble() %>% 
  gather(key = "key", value = "value", -var)


pgpea <- ggplot()+
  geom_vline(xintercept = 0, col="grey70") +
  stat_density_ridges(data  = posteriorsimpelmodellong, 
                      aes(x      = value,
                          y      = var,
                          fill   = as.factor(var),
                          vline_color = ..quantile..),
                      scale = 1.8, 
                      alpha = .8,
                      size=0.5,
                      rel_min_height=.01,
                      quantile_lines = TRUE, 
                      quantiles = c(0.025,0.5,0.975)) +
  scale_discrete_manual("vline_color",
                        values = c("blue", "red", "blue", "black"), 
                        breaks = c(1, 2),
                        labels = c("5% & 95% quantiles", "mean"),
                        name = NULL) +
  scale_y_discrete(labels=c("beta_PEA_Caceres"=parse(text = TeX(paste0("$\\beta_{",pea,",Caceres}$"))),
                            "beta_PEA_Madrid"=parse(text = TeX(paste0("$\\beta_{",pea,",Madrid}$"))),
                            "beta_PEA_Portugal"=parse(text = TeX(paste0("$\\beta_{",pea,",Portugal}$"))),
                            "beta_PEA_Asturias"=parse(text = TeX(paste0("$\\beta_{",pea,",Asturias}$"))),
                            "beta_PEA_Bordeaux"=parse(text = TeX(paste0("$\\beta_{",pea,",Bordeaux}$")))
  )) + 
  labs(fill = "Sites",
       y = "", 
       x="") +
  scale_fill_manual(values=c(vir_lite("cyan2",ds=ds),
                             vir_lite("navyblue",ds=ds),
                             vir_lite("pink",ds=ds),
                             vir_lite("deeppink",ds=ds),
                             vir_lite("dodgerblue2",ds=ds)),
                    labels = c("Asturias", 
                               "Bordeaux",
                               "Caceres",
                               "Madrid",
                               "Portugal")) +
  theme_bw() + theme(axis.text = element_text(size=22),axis.title = element_text(size=22), 
                     legend.text = element_text(size=16),legend.title = element_text(size=17),
                     plot.title = element_text(size=22))   + guides(vline_color = FALSE)



# >> M12 => rPEAs    ####
# ---------------    "

mod <- readRDS(file=paste0("outputs/models/",part,"/MOD12.rds"))
pea="rPEA"
variable ="rPEA.sc"

POST <- posterior_samples(mod,pars = "^r_site\\[.*sc\\]") %>% dplyr::rename(
  beta_PEA_Portugal = paste0('r_site[portugal,',variable,']'),
  beta_PEA_Bordeaux = paste0('r_site[bordeaux,',variable,']'),
  beta_PEA_Asturias = paste0('r_site[asturias,',variable,']'), 
  beta_PEA_Madrid = paste0('r_site[madrid,',variable,']'),
  beta_PEA_Caceres = paste0('r_site[caceres,',variable,']')
)
POST <- as.data.frame(t(POST))
POST$var <- as.factor(rownames(POST))

posteriorsimpelmodellong <- POST %>%  as_tibble() %>% 
  gather(key = "key", value = "value", -var)


prpea <- ggplot()+
  geom_vline(xintercept = 0, col="grey70") +
  stat_density_ridges(data  = posteriorsimpelmodellong, 
                      aes(x      = value,
                          y      = var,
                          fill   = as.factor(var),
                          vline_color = ..quantile..),
                      scale = 1.8, 
                      alpha = .8,
                      size=0.5,
                      rel_min_height=.01,
                      quantile_lines = TRUE, 
                      quantiles = c(0.025,0.5,0.975)) +
  scale_discrete_manual("vline_color",
                        values = c("blue", "red", "blue", "black"), 
                        breaks = c(1, 2),
                        labels = c("5% & 95% quantiles", "mean"),
                        name = NULL) +
  scale_y_discrete(labels=c("beta_PEA_Caceres"=parse(text = TeX(paste0("$\\beta_{",pea,",Caceres}$"))),
                            "beta_PEA_Madrid"=parse(text = TeX(paste0("$\\beta_{",pea,",Madrid}$"))),
                            "beta_PEA_Portugal"=parse(text = TeX(paste0("$\\beta_{",pea,",Portugal}$"))),
                            "beta_PEA_Asturias"=parse(text = TeX(paste0("$\\beta_{",pea,",Asturias}$"))),
                            "beta_PEA_Bordeaux"=parse(text = TeX(paste0("$\\beta_{",pea,",Bordeaux}$")))
  )) + 
  labs(title="",
       y        = "", 
       x=""
  ) + 
  scale_fill_manual(values=c(vir_lite("cyan2",ds=ds),
                             vir_lite("navyblue",ds=ds),
                             vir_lite("pink",ds=ds),
                             vir_lite("deeppink",ds=ds),
                             vir_lite("dodgerblue2",ds=ds))) +
  theme_bw() + theme(axis.text = element_text(size=22),axis.title = element_text(size=22), 
                     legend.position = "none", plot.title = element_text(size=22)) 
p1 <- ggarrange(pGP,pCP,labels=c("M9","M10"),font.label = list(size = 20),nrow=1,widths = c(1.2,1))
p2 <- ggarrange(pgpea,prpea,labels=c("M11","M12"),font.label = list(size = 20),nrow=1,widths = c(1.2,1))
fig <- ggarrange(p1,p2,nrow=2)



# >> Figures S16, S18 and S20 ####

if(part=="P1"){
  ggsave(fig, file=paste0("figs/SuppInfo/",part,"M9toM12PosteriorDistri.png"),width=20,height=12)  
} else{
  ggsave(fig, file=paste0("figs/SuppInfo/",part,"M9toM12PosteriorDistri.png"),width=20,height=12)  
}


