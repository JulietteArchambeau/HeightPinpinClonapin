################################################################################################"
###########      Script to create the figure S19 of the Supplementary Information      #########"
################################################################################################"


# This figure summarize the genetic response in models M3, M5 & M6


library(cowplot)
library(ggplot2)


data <- readRDS(file="data/TrainP1.RDS")




# Model M3 ####

# The genetic response is represented by the genotype and provenance intercepts.

# a) Provenace intercepts ####

mod <- readRDS(file="outputs/models/P1/MOD3.rds")
POST <- posterior_samples(mod,pars = "^r_prov\\[")
colnames(POST) <- str_sub(colnames(POST),8,-12)
POST <- as.data.frame(t(POST))
POST$prov <- as.factor(rownames(POST))

data <- droplevels(data)
ps <- data %>%
  group_by(prov) %>% 
  summarise_at(vars(paste0(rep("Q",6),1:6)), mean)
ps$max.Q.prov <- colnames(ps[,2:7])[apply(ps[,2:7],1,which.max)]
ps$prov <- as.factor(ps$prov)

posteriorsimpelmodellong <- inner_join(POST, ps[,c("prov","max.Q.prov")],by="prov") %>%  as_tibble() %>% 
  gather(key = "key", value = "value", -prov,-max.Q.prov)%>%
  group_by(prov) %>%
  dplyr::mutate(meanperprov = mean(value))%>%
  ungroup()

pM3_prov <- ggplot()+
  stat_density_ridges(data  = posteriorsimpelmodellong, 
                      aes(x      = value,
                          y      = reorder(as.factor(prov), meanperprov),
                          height = ..density.., 
                          fill   = as.factor(max.Q.prov),
                          vline_color = ..quantile..),
                      scale = 3, 
                      alpha = .6,
                      rel_min_height=c(.01),
                      size=0.2,
                      quantile_lines = TRUE, quantiles = c(0.025,0.5,0.975)) +
  
  scale_y_discrete(labels=c("ALT"=parse(text = TeX("$ALT$")),
                            "ARM"=parse(text = TeX("$ARM$")),
                            "ARN"=parse(text = TeX("$ARN$")),
                            "BAY"=parse(text = TeX("$BAY$")),
                            "BON"=parse(text = TeX("$BON$")),
                            "CAD"=parse(text = TeX("$CAD$")),
                            "CAR"=parse(text = TeX("$CAR$")),
                            "CAS"=parse(text = TeX("$CAS$")),
                            "CEN"=parse(text = TeX("$CEN$")),
                            "COC"=parse(text = TeX("$COC$")),
                            "COM"=parse(text = TeX("$COM$")),
                            "CUE"=parse(text = TeX("$CUE$")),
                            "HOU"=parse(text = TeX("$HOU$")),
                            "LAM"=parse(text = TeX("$LAM$")),
                            "LEI"=parse(text = TeX("$LEI$")),
                            "MAD"=parse(text = TeX("$MAD$")),
                            "MIM"=parse(text = TeX("$MIM$")),
                            "OLB"=parse(text = TeX("$OLB$")),
                            "OLO"=parse(text = TeX("$OLO$")),
                            "ORI"=parse(text = TeX("$ORI$")),
                            "PET"=parse(text = TeX("$PET$")),
                            "PIA"=parse(text = TeX("$PIA$")),
                            "PIE"=parse(text = TeX("$PIE$")),
                            "PLE"=parse(text = TeX("$PLE$")),
                            "PUE"=parse(text = TeX("$PUE$")),
                            "QUA"=parse(text = TeX("$QUA$")),
                            "SAC"=parse(text = TeX("$SAC$")),
                            "SAL"=parse(text = TeX("$SAL$")),
                            "SEG"=parse(text = TeX("$SEG$")),
                            "SIE"=parse(text = TeX("$SIE$")),
                            "STJ"=parse(text = TeX("$STJ$")),
                            "TAM"=parse(text = TeX("$TAM$")),
                            "VAL"=parse(text = TeX("$VAL$")),
                            "VER"=parse(text = TeX("$VER$")))) +
  
  coord_cartesian(c(-0.35,0.25),c(0.5,37))+
  
  scale_discrete_manual("vline_color",
                        values = c("blue", "red", "blue", "black"), 
                        breaks = c(1, 2),
                        labels = c("2.5 and 97.5th percentiles", "Median"),
                        name = NULL) +
  scale_fill_manual(values=c("orangered3",
                             "gold2","darkorchid3",
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
       x=""
       #x        = TeX("Provenance intercept P_{p}")
  )+
  theme_bw() + theme(axis.text = element_text(size=12),
                     axis.title = element_text(size=14), 
                     legend.text = element_text(size=18),
                     legend.title = element_text(size=20))

#---------------------------------------------------------------------------------------------------"


# Model M4 ####

# The genetic response is represented by the genotype and provenance intercepts, and the gene pool intercepts.

mod <- readRDS(file="outputs/models/P1/MOD4.rds")

## a) Provenance intercepts ####

POST <- posterior_samples(mod,pars = "^r_prov\\[")
colnames(POST) <- str_sub(colnames(POST),8,-12)
POST <- as.data.frame(t(POST))
POST$prov <- as.factor(rownames(POST))

data <- droplevels(data)
ps <- data %>%
  group_by(prov) %>% 
  summarise_at(vars(paste0(rep("Q",6),1:6)), mean)
ps$max.Q.prov <- colnames(ps[,2:7])[apply(ps[,2:7],1,which.max)]
ps$prov <- as.factor(ps$prov)

posteriorsimpelmodellong <- inner_join(POST, ps[,c("prov","max.Q.prov")],by="prov") %>%  as_tibble() %>% 
  gather(key = "key", value = "value", -prov,-max.Q.prov)%>%
  group_by(prov) %>%
  dplyr::mutate(meanperprov = mean(value))%>%
  ungroup()

pM4_prov <- ggplot() +
  stat_density_ridges(data  = posteriorsimpelmodellong, 
                      aes(x      = value,
                          y      = reorder(as.factor(prov), meanperprov),
                          height = ..density.., 
                          fill   = as.factor(max.Q.prov),
                          vline_color = ..quantile..),
                      scale = 2, 
                      alpha = .6,
                      rel_min_height=c(.01),
                      size=0.2,
                      quantile_lines = TRUE, quantiles = c(0.025,0.5,0.975)) +
  
  scale_y_discrete(labels=c("ALT"=parse(text = TeX("$ALT$")),
                            "ARM"=parse(text = TeX("$ARM$")),
                            "ARN"=parse(text = TeX("$ARN$")),
                            "BAY"=parse(text = TeX("$BAY$")),
                            "BON"=parse(text = TeX("$BON$")),
                            "CAD"=parse(text = TeX("$CAD$")),
                            "CAR"=parse(text = TeX("$CAR$")),
                            "CAS"=parse(text = TeX("$CAS$")),
                            "CEN"=parse(text = TeX("$CEN$")),
                            "COC"=parse(text = TeX("$COC$")),
                            "COM"=parse(text = TeX("$COM$")),
                            "CUE"=parse(text = TeX("$CUE$")),
                            "HOU"=parse(text = TeX("$HOU$")),
                            "LAM"=parse(text = TeX("$LAM$")),
                            "LEI"=parse(text = TeX("$LEI$")),
                            "MAD"=parse(text = TeX("$MAD$")),
                            "MIM"=parse(text = TeX("$MIM$")),
                            "OLB"=parse(text = TeX("$OLB$")),
                            "OLO"=parse(text = TeX("$OLO$")),
                            "ORI"=parse(text = TeX("$ORI$")),
                            "PET"=parse(text = TeX("$PET$")),
                            "PIA"=parse(text = TeX("$PIA$")),
                            "PIE"=parse(text = TeX("$PIE$")),
                            "PLE"=parse(text = TeX("$PLE$")),
                            "PUE"=parse(text = TeX("$PUE$")),
                            "QUA"=parse(text = TeX("$QUA$")),
                            "SAC"=parse(text = TeX("$SAC$")),
                            "SAL"=parse(text = TeX("$SAL$")),
                            "SEG"=parse(text = TeX("$SEG$")),
                            "SIE"=parse(text = TeX("$SIE$")),
                            "STJ"=parse(text = TeX("$STJ$")),
                            "TAM"=parse(text = TeX("$TAM$")),
                            "VAL"=parse(text = TeX("$VAL$")),
                            "VER"=parse(text = TeX("$VER$")))) +
  
  coord_cartesian(c(-0.35,0.25),c(0.5,36))+
  
  scale_discrete_manual("vline_color",
                        values = c("blue", "red", "blue", "black"), 
                        breaks = c(1, 2),
                        labels = c("2.5 and 97.5th percentiles", "Median"),
                        name = NULL) +
  
  scale_fill_manual(values=c("orangered3",
                             "gold2",
                             "darkorchid3",
                             "navyblue",
                             "turquoise2",
                             "green3"), labels = c("Northern Africa", 
                                                   "Corsica",
                                                   "Central Spain",
                                                   "French Atlantic",
                                                   "Iberian Atlantic",
                                                   "South-eastern Spain"),name="Gene pools") +
  labs(fill     = "Gene pools",
       y        = "", 
       x=""
       #x        = TeX("Provenance intercept P_{p}")
  ) + 
  theme_bw() + 
  theme(axis.text = element_text(size=12),axis.title = element_text(size=14), 
        legend.text = element_text(size=18),legend.title = element_text(size=20))


# b) Gene pool intercepts ####
POST <- posterior_samples(mod,pars = "^r_mmQ1Q2Q3Q4Q5Q6\\[")
colnames(POST) <- str_sub(colnames(POST),18,-12)
POST <- as.data.frame(t(POST))
POST$genepool <- as.factor(rownames(POST))

posteriorsimpelmodellong <- POST %>%  as_tibble() %>% 
  gather(key = "key", value = "value", -genepool)%>%
  group_by(genepool) %>%
  dplyr::mutate(meanpergenepool = mean(value))%>%
  ungroup()

pM4_GP <- ggplot()+
  geom_vline(xintercept = 0, 
             col        = "grey70") +
  stat_density_ridges(data  = posteriorsimpelmodellong, 
                      aes(x      = value,
                          y      = reorder(as.factor(genepool), meanpergenepool),
                          # height = ..density.., 
                          fill   = as.factor(genepool),
                          vline_color = ..quantile..),
                      scale = 2, 
                      alpha = .6,
                      rel_min_height=c(.001),
                      size=0.2,
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
                             "green3"), labels = c("Northern Africa", 
                                                   "Corsica",
                                                   "Central Spain",
                                                   "French Atlantic",
                                                   "Iberian Atlantic",
                                                   "South-eastern Spain")) +
  
  
  scale_y_discrete(labels=c("Q1"=parse(text = TeX("$NA$")),
                            "Q2"=parse(text = TeX("$C$")),
                            "Q3"=parse(text = TeX("$CS$")),
                            "Q4"=parse(text = TeX("$FA$")),
                            "Q5"=parse(text = TeX("$IA$")),
                            "Q6"=parse(text = TeX("$SES$")))) + 
  labs(fill     = "Gene pools",
       y        = "", 
       x=""
       #x        = TeX("Gene pool intercept g_{j}")
  ) + 
  theme_bw() + theme(axis.text = element_text(size=12),axis.title = element_text(size=14), 
                     legend.text = element_text(size=18),legend.title = element_text(size=20))


#---------------------------------------------------------------------------------------------------"

# Model M6 ####

# The genetic response is represented by the genotype and provenance intercepts, the gene pool intercepts and the provenance climate-of-origin intercepts. 

mod <- readRDS(file="outputs/models/P1/MOD6.rds")

# a) Provenance intercepts ####
POST <- posterior_samples(mod,pars = "^r_prov\\[")
colnames(POST) <- str_sub(colnames(POST),8,-12)
POST <- as.data.frame(t(POST))
POST$prov <- as.factor(rownames(POST))

data <- droplevels(data)
ps <- data %>%
  group_by(prov) %>% 
  summarise_at(vars(paste0(rep("Q",6),1:6)), mean)
ps$max.Q.prov <- colnames(ps[,2:7])[apply(ps[,2:7],1,which.max)]
ps$prov <- as.factor(ps$prov)

posteriorsimpelmodellong <- inner_join(POST, ps[,c("prov","max.Q.prov")],by="prov") %>%  as_tibble() %>% 
  gather(key = "key", value = "value", -prov,-max.Q.prov)%>%
  group_by(prov) %>%
  dplyr::mutate(meanperprov = mean(value))%>%
  ungroup()

pM6_prov <- ggplot()+
  stat_density_ridges(data  = posteriorsimpelmodellong, 
                      aes(x      = value,
                          y      = reorder(as.factor(prov), meanperprov),
                          height = ..density.., 
                          fill   = as.factor(max.Q.prov),
                          vline_color = ..quantile..),
                      scale = 3, 
                      alpha = .6,
                      rel_min_height=c(.01),
                      size=0.2,
                      quantile_lines = TRUE, quantiles = c(0.025,0.5,0.975)) +
  
  coord_cartesian(c(-0.35,0.25),c(0.5,36))+
  
  scale_discrete_manual("vline_color",
                        values = c("blue", "red", "blue", "black"), 
                        breaks = c(1, 2),
                        labels = c("2.5 and 97.5th percentiles", "Median"),
                        name = NULL) +
  
  
  scale_fill_manual(values=c("orangered3",
                             "gold2",
                             "darkorchid3",
                             "navyblue",
                             "turquoise2",
                             "green3"), labels = c("Northern Africa", 
                                                   "Corsica",
                                                   "Central Spain",
                                                   "French Atlantic",
                                                   "Iberian Atlantic",
                                                   "South-eastern Spain")) +
  
  scale_y_discrete(labels=c("ALT"=parse(text = TeX("$ALT$")),
                            "ARM"=parse(text = TeX("$ARM$")),
                            "ARN"=parse(text = TeX("$ARN$")),
                            "BAY"=parse(text = TeX("$BAY$")),
                            "BON"=parse(text = TeX("$BON$")),
                            "CAD"=parse(text = TeX("$CAD$")),
                            "CAR"=parse(text = TeX("$CAR$")),
                            "CAS"=parse(text = TeX("$CAS$")),
                            "CEN"=parse(text = TeX("$CEN$")),
                            "COC"=parse(text = TeX("$COC$")),
                            "COM"=parse(text = TeX("$COM$")),
                            "CUE"=parse(text = TeX("$CUE$")),
                            "HOU"=parse(text = TeX("$HOU$")),
                            "LAM"=parse(text = TeX("$LAM$")),
                            "LEI"=parse(text = TeX("$LEI$")),
                            "MAD"=parse(text = TeX("$MAD$")),
                            "MIM"=parse(text = TeX("$MIM$")),
                            "OLB"=parse(text = TeX("$OLB$")),
                            "OLO"=parse(text = TeX("$OLO$")),
                            "ORI"=parse(text = TeX("$ORI$")),
                            "PET"=parse(text = TeX("$PET$")),
                            "PIA"=parse(text = TeX("$PIA$")),
                            "PIE"=parse(text = TeX("$PIE$")),
                            "PLE"=parse(text = TeX("$PLE$")),
                            "PUE"=parse(text = TeX("$PUE$")),
                            "QUA"=parse(text = TeX("$QUA$")),
                            "SAC"=parse(text = TeX("$SAC$")),
                            "SAL"=parse(text = TeX("$SAL$")),
                            "SEG"=parse(text = TeX("$SEG$")),
                            "SIE"=parse(text = TeX("$SIE$")),
                            "STJ"=parse(text = TeX("$STJ$")),
                            "TAM"=parse(text = TeX("$TAM$")),
                            "VAL"=parse(text = TeX("$VAL$")),
                            "VER"=parse(text = TeX("$VER$")))) +
  labs(fill     = "Gene pools",
       y        = "", 
       x        = TeX("Provenance intercept P_{p}")
  ) + 
  theme_bw() + theme(axis.text = element_text(size=12),
                     axis.title = element_text(size=16), 
                     legend.text = element_text(size=18),
                     legend.title = element_text(size=20))


## b) Gene pool intercepts ####

POST <- posterior_samples(mod,pars = "^r_mmQ1Q2Q3Q4Q5Q6\\[")
colnames(POST) <- str_sub(colnames(POST),18,-12)
POST <- as.data.frame(t(POST))
POST$genepool <- as.factor(rownames(POST))

posteriorsimpelmodellong <- POST %>%  as_tibble() %>% 
  gather(key = "key", value = "value", -genepool)%>%
  group_by(genepool) %>%
  dplyr::mutate(meanpergenepool = mean(value))%>%
  ungroup()

pM6_GP <- ggplot()+
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
                      size=0.2,
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
                             "green3"), labels = c("Northern Africa", 
                                                   "Corsica",
                                                   "Central Spain",
                                                   "French Atlantic",
                                                   "Iberian Atlantic",
                                                   "South-eastern Spain")) +
  
  
  scale_y_discrete(labels=c("Q1"=parse(text = TeX("$NA$")),
                            "Q2"=parse(text = TeX("$C$")),
                            "Q3"=parse(text = TeX("$CS$")),
                            "Q4"=parse(text = TeX("$FA$")),
                            "Q5"=parse(text = TeX("$IA$")),
                            "Q6"=parse(text = TeX("$SES$")))) + 
  labs(fill     = "Gene pools",
       y        = "", 
       x        = TeX("Gene pool intercept g_{j}")
  ) + 
  theme_bw() + theme(axis.text = element_text(size=12),
                     axis.title = element_text(size=16), 
                     legend.text = element_text(size=18),
                     legend.title = element_text(size=20))


## c) Provenance climate intercepts ####

POST <- posterior_samples(mod,pars = "^r_prov_clim")
colnames(POST) <- str_sub(colnames(POST),13,-12)
POST <- as.data.frame(t(POST))
POST$prov <- as.factor(rownames(POST))

data <- droplevels(data)
ps <- data %>%
  group_by(prov) %>% 
  summarise_at(vars(paste0(rep("Q",6),1:6)), mean)
ps$max.Q.prov <- colnames(ps[,2:7])[apply(ps[,2:7],1,which.max)]
ps$prov <- as.factor(ps$prov)

posteriorsimpelmodellong <- inner_join(POST, ps[,c("prov","max.Q.prov")],by="prov") %>%  as_tibble() %>% 
  gather(key = "key", value = "value", -prov,-max.Q.prov)%>%
  group_by(prov) %>%
  dplyr::mutate(meanperprov = mean(value))%>%
  ungroup()


pM6_PC <- ggplot()+
  stat_density_ridges(data  = posteriorsimpelmodellong, 
                      aes(x      = value,
                          y      = reorder(as.factor(prov), meanperprov),
                          height = ..density.., 
                          fill   = as.factor(max.Q.prov),
                          vline_color = ..quantile..),
                      scale = 3, 
                      alpha = .6,
                      rel_min_height=c(.01),
                      size=0.2,
                      quantile_lines = TRUE, quantiles = c(0.025,0.5,0.975)) +
  
  
  scale_y_discrete(labels=c("ALT"=parse(text = TeX("$ALT$")),
                            "ARM"=parse(text = TeX("$ARM$")),
                            "ARN"=parse(text = TeX("$ARN$")),
                            "BAY"=parse(text = TeX("$BAY$")),
                            "BON"=parse(text = TeX("$BON$")),
                            "CAD"=parse(text = TeX("$CAD$")),
                            "CAR"=parse(text = TeX("$CAR$")),
                            "CAS"=parse(text = TeX("$CAS$")),
                            "CEN"=parse(text = TeX("$CEN$")),
                            "COC"=parse(text = TeX("$COC$")),
                            "COM"=parse(text = TeX("$COM$")),
                            "CUE"=parse(text = TeX("$CUE$")),
                            "HOU"=parse(text = TeX("$HOU$")),
                            "LAM"=parse(text = TeX("$LAM$")),
                            "LEI"=parse(text = TeX("$LEI$")),
                            "MAD"=parse(text = TeX("$MAD$")),
                            "MIM"=parse(text = TeX("$MIM$")),
                            "OLB"=parse(text = TeX("$OLB$")),
                            "OLO"=parse(text = TeX("$OLO$")),
                            "ORI"=parse(text = TeX("$ORI$")),
                            "PET"=parse(text = TeX("$PET$")),
                            "PIA"=parse(text = TeX("$PIA$")),
                            "PIE"=parse(text = TeX("$PIE$")),
                            "PLE"=parse(text = TeX("$PLE$")),
                            "PUE"=parse(text = TeX("$PUE$")),
                            "QUA"=parse(text = TeX("$QUA$")),
                            "SAC"=parse(text = TeX("$SAC$")),
                            "SAL"=parse(text = TeX("$SAL$")),
                            "SEG"=parse(text = TeX("$SEG$")),
                            "SIE"=parse(text = TeX("$SIE$")),
                            "STJ"=parse(text = TeX("$STJ$")),
                            "TAM"=parse(text = TeX("$TAM$")),
                            "VAL"=parse(text = TeX("$VAL$")),
                            "VER"=parse(text = TeX("$VER$")))) +
  
  coord_cartesian(c(-0.35,0.25),c(0.5,35))+
  
  scale_discrete_manual("vline_color",
                        values = c("blue", "red", "blue", "black"), 
                        breaks = c(1, 2),
                        labels = c("2.5 and 97.5th percentiles", "Median"),
                        name = NULL) +
  
  scale_fill_manual(values=c("orangered3",
                             "gold2",
                             "darkorchid3",
                             "navyblue",
                             "turquoise2",
                             "green3"), labels = c("Northern Africa", 
                                                   "Corsica",
                                                   "Central Spain",
                                                   "French Atlantic",
                                                   "Iberian Atlantic",
                                                   "South-eastern Spain"),name="Gene pools") +
  labs(fill     = "Gene pools",
       y        = "", 
       x        = TeX("Provenance climate intercepts cp_{p}")
       ) + 
  theme_bw() + theme(axis.text = element_text(size=12),
                     axis.title = element_text(size=16), 
                     legend.text = element_text(size=18),
                     legend.title = element_text(size=20))



## Final figure ####


# Extract legend
legend <- get_legend(pM3_prov)


# Then remove legend
pM3_prov <- pM3_prov + theme(legend.position = "none")

pM4_prov <- pM4_prov + theme(legend.position = "none") 
pM4_GP <- pM4_GP + theme(legend.position = "none") 

pM6_prov <- pM6_prov + theme(legend.position = "none")
pM6_GP <- pM6_GP + theme(legend.position = "none")
pM6_PC <- pM6_PC + theme(legend.position = "none")



fig <- plot_grid(pM3_prov,NULL,legend,pM4_prov,pM4_GP,NULL,pM6_prov,pM6_GP,pM6_PC, labels = c('M3','',"","M4","","","M6","",""),
                 label_size = 20 )
ggsave(fig, file="figs/SuppInfo/PostDHPC.png",width = 16,height=19) 
