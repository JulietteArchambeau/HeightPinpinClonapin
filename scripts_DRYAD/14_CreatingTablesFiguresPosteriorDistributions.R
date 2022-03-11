########################################################################################################################"
#                                                                                                                      #
#               Tables and figures of the posterior distributions in the Supplementary Information                     #
#                                                                                                                      #
#                                   Juliette Archambeau                                                                #
#                                       11/03/2022                                                                     #
#                                                                                                                      #
########################################################################################################################"

library(broom)
library(latex2exp)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(stringi)
library(tidybayes)
library(dplyr)
library(bayesplot)
color_scheme_set("green")
library(devtools)
library(xtable)
library(ggridges)
library(tidyverse)
library(tibble)
library(brms)


# Functions used in the script:
source("scripts/Functions/vir_lite.R") # available here: https://github.com/JulietteArchambeau/HeightPinpinClonapin/blob/master/scripts/Functions/vir_lite.R
square <- function(x) (x*x)
ds=0.7 # Parameter of the funtion vir_lite


# In the Supplementary Information, we report the medians and the 95% credible intervals of the posterior distributions.
prob=0.95
probs <- c((1 - prob) / 2, 1 - (1 - prob) / 2)

# Load the train dataset of the P1 partition
#data <- readRDS(file="../../data/TrainP1.RDS")


# Data partition
part <- "P1" # choose between P1, P2 and P3


# Model M0 ####
# ======== "

mod <- readRDS(file= paste0("outputs/models/",part,"/MOD0.rds"))


# >> Tables S14 and S36. ####
# ---------------------- "


# extract the standard deviations:
pars <- mod %>% get_variables() %>%  as.vector() %>% str_subset("^(sd|sigma)")
sims <- as.array(mod, pars = pars, fixed = TRUE)

df <- as.data.frame(matrix(NA,length(pars),4,dimnames = list(pars,c("Median","SD","InfCI","SupCI"))))
for(i in pars){
  sims_i <- square(sims[, , i]) # we want the variances and not the standard deviations
  quan <- unname(quantile(sims_i, probs = probs))
  df[i,"InfCI"] <- quan[1]
  df[i,"SupCI"] <- quan[2]
  df[i,"Median"] <- median(sims_i)
  df[i,"SD"] <- sd(sims_i)}

df <- df %>% mutate(Parameter = recode_factor(rownames(df),'sigma'='$\\sigma^{2}$',
                                              'sd_prov__Intercept'="$\\sigma^{2}_{P}$",
                                              'sd_prov:clon__Intercept'='$\\sigma^{2}_{G}$',
                                              'sd_site__Intercept'='$\\sigma^{2}_{S}$',
                                              'sd_site:block__Intercept'='$\\sigma^{2}_{B}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)

# Extract the coefficients of the fixed effects
pars <- mod %>% get_variables() %>%  as.vector() %>% str_subset("^(b)")
sims <- as.array(mod, pars = pars, fixed = TRUE)

df2 <- as.data.frame(matrix(NA,length(pars),4,dimnames = list(pars,c("Median","SD","InfCI","SupCI"))))
for(i in pars){
  sims_i <- sims[, , i]
  quan <- unname(quantile(sims_i, probs = probs))
  df2[i,"InfCI"] <- quan[1]
  df2[i,"SupCI"] <- quan[2]
  df2[i,"Median"] <- median(sims_i)
  df2[i,"SD"] <- sd(sims_i)}

df2 <- df2 %>% 
  mutate(Parameter = recode_factor(rownames(df2),'b_age.sc'="$\\beta_{age}$",
                                                'b_Iage.scE2'= '$\\beta_{age2}$',
                                                'b_Intercept'='$\\beta_{0}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)


tab <- bind_rows(df,df2)

print(xtable(tab, type = "latex",digits=3), 
      file = paste0("tables/Posteriors/M0_MainVarPost.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})


# Model M1 ####
# ======== "

mod <- readRDS(file=paste0("outputs/models/",part,"/MOD1.rds"))



# >> Tables S15 and S37. ####
# ---------------------- "

pars <- mod %>% get_variables() %>%  as.vector() %>% str_subset("^(sd|sigma)")
sims <- as.array(mod, pars = pars, fixed = TRUE)

df <- as.data.frame(matrix(NA,length(pars),4,dimnames = list(pars,c("Median","SD","InfCI","SupCI"))))
for(i in pars){
  sims_i <- square(sims[, , i])
  quan <- unname(quantile(sims_i, probs = probs))
  df[i,"InfCI"] <- quan[1]
  df[i,"SupCI"] <- quan[2]
  df[i,"Median"] <- median(sims_i)
  df[i,"SD"] <- sd(sims_i)}

df <- df %>% mutate(Parameter = recode_factor(rownames(df),
                                              'sigma'='$\\sigma^{2}$',
                                              'sd_prov__Intercept'="$\\sigma^{2}_{P}$",
                                              'sd_prov:clon__Intercept'='$\\sigma^{2}_{G}$',
                                              'sd_site__Intercept'='$\\sigma^{2}_{S}$',
                                              'sd_site:block__Intercept'='$\\sigma^{2}_{B}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)


pars <- mod %>% get_variables() %>%  as.vector() %>% str_subset("^(b)")
sims <- as.array(mod, pars = pars, fixed = TRUE)

df2 <- as.data.frame(matrix(NA,length(pars),4,dimnames = list(pars,c("Median","SD","InfCI","SupCI"))))
for(i in pars){
  sims_i <- sims[, , i]
  quan <- unname(quantile(sims_i, probs = probs))
  df2[i,"InfCI"] <- quan[1]
  df2[i,"SupCI"] <- quan[2]
  df2[i,"Median"] <- median(sims_i)
  df2[i,"SD"] <- sd(sims_i)}
df2 <- df2 %>% mutate(Parameter = recode_factor(rownames(df2),
                                                'b_age.sc'="$\\beta_{age}$",
                                                'b_Iage.scE2'= '$\\beta_{age2}$',
                                                'b_Intercept'='$\\beta_{0}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)


tab <- bind_rows(df,df2)
print(xtable(tab, type = "latex",digits=3), 
      file = paste0("tables/Posteriors/M1_MainVarPost.tex"), include.rownames=FALSE,sanitize.text.function = function(x) {x})


# >> Table S16. ####
# ------------- "

# only for the P1 partition

df <- mod %>% broom::tidyMCMC(estimate.method = "median",conf.int = T,conf.level = 0.95) %>% 
  filter(str_detect(term, "^(r_site\\[)")) %>% 
  rename_all(str_to_title) %>% 
  dplyr::rename("Median"=Estimate,"SD"=Std.error,"InfCI"=Conf.low,"SupCI"=Conf.high) %>% 
  mutate(Parameter = recode_factor(Term,
                                   'r_site[asturias,Intercept]'="$S_{Asturias}$",
                                   'r_site[bordeaux,Intercept]'= '$S_{Bordeaux}$',
                                   'r_site[caceres,Intercept]'='$S_{Caceres}$',
                                   'r_site[madrid,Intercept]'='$S_{Madrid}$',
                                   'r_site[portugal,Intercept]'="$S_{Portugal}$")) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)

print(xtable(df, type = "latex",digits=3), 
      file = paste0("tables/Posteriors/M1_SiteInterceptsPost.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})



# >> Table S17. ####
# ------------- "

# only for the P1 partition

df <- mod %>% broom::tidyMCMC(estimate.method = "mean",conf.int = T) %>% 
  filter(str_detect(term, "^(r_prov\\[)")) %>% 
  rename_all(str_to_title) %>% 
  dplyr::rename("Median"=Estimate,"SD"=Std.error,"InfCI"=Conf.low,"SupCI"=Conf.high) %>% 
  mutate(Parameter= recode_factor(Term, 
                                  'r_prov[MIM,Intercept]' = '$P_{MIM}$',
                                  'r_prov[CEN,Intercept]' = '$P_{CEN}$',
                                  'r_prov[ORI,Intercept]' = '$P_{ORI}$',
                                  'r_prov[STJ,Intercept]' = '$P_{STJ}$',
                                  'r_prov[HOU,Intercept]' = '$P_{HOU}$',
                                  'r_prov[CUE,Intercept]' = '$P_{CUE}$',
                                  'r_prov[TAM,Intercept]' = '$P_{TAM}$',
                                  'r_prov[LEI,Intercept]' = '$P_{LEI}$',
                                  'r_prov[VER,Intercept]' = '$P_{VER}$',
                                  'r_prov[CAS,Intercept]' = '$P_{CAS}$',
                                  'r_prov[PIA,Intercept]' = '$P_{PIA}$',
                                  'r_prov[ARM,Intercept]' = '$P_{ARM}$',
                                  'r_prov[PET,Intercept]' = '$P_{PET}$',
                                  'r_prov[VAL,Intercept]' = '$P_{VAL}$',
                                  'r_prov[SAL,Intercept]' = '$P_{SAL}$',
                                  'r_prov[OLO,Intercept]' = '$P_{OLO}$',
                                  'r_prov[CAD,Intercept]' = '$P_{CAD}$',
                                  'r_prov[ARN,Intercept]' = '$P_{ARN}$',
                                  'r_prov[BAY,Intercept]' = '$P_{BAY}$',
                                  'r_prov[SIE,Intercept]' = '$P_{SIE}$',
                                  'r_prov[SEG,Intercept]' = '$P_{SEG}$',
                                  'r_prov[PLE,Intercept]' = '$P_{PLE}$',
                                  'r_prov[BON,Intercept]' = '$P_{BON}$',
                                  'r_prov[COC,Intercept]' = '$P_{COC}$',
                                  'r_prov[SAC,Intercept]' = '$P_{SAC}$',
                                  'r_prov[QUA,Intercept]' = '$P_{QUA}$',
                                  'r_prov[CAR,Intercept]' = '$P_{CAR}$',
                                  'r_prov[OLB,Intercept]' = '$P_{OLB}$',
                                  'r_prov[PIE,Intercept]' = '$P_{PIE}$',
                                  'r_prov[PUE,Intercept]' = '$P_{PUE}$',
                                  'r_prov[ALT,Intercept]' = '$P_{ALT}$',
                                  'r_prov[LAM,Intercept]' = '$P_{LAM}$',
                                  'r_prov[MAD,Intercept]' = '$P_{MAD}$',
                                  'r_prov[COM,Intercept]' = '$P_{COM}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)


print(xtable(df, type = "latex",digits=3), 
      file = paste0("tables/Posteriors/M1_ProvInterceptsPost.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})


# >> Figure S11. ####
# -------------  "

# only for the P1 partition

p <- plot(conditional_effects(mod,"age.sc"),plot=FALSE)[[1]] + 
  xlab("Mean-centered age") + 
  ylab("Logarithm of height (mm)") +
  theme_bw()

ggsave(p,file="figs/SuppInfo/M1_CondEffectAge.png",height=6,width=6)



# Model M2 ####
# ======== "

mod <- readRDS(file= paste0("outputs/models/",part,"/MOD2.rds"))


# >> Tables S18 and S38. ####
# ---------------------- "


pars <- mod %>% get_variables() %>%  as.vector() %>% str_subset("^(sd|sigma)")
sims <- as.array(mod, pars = pars, fixed = TRUE)

df <- as.data.frame(matrix(NA,length(pars),4,dimnames = list(pars,c("Median","SD","InfCI","SupCI"))))
for(i in pars){
  sims_i <- square(sims[, , i])
  quan <- unname(quantile(sims_i, probs = probs))
  df[i,"InfCI"] <- quan[1]
  df[i,"SupCI"] <- quan[2]
  df[i,"Median"] <- median(sims_i)
  df[i,"SD"] <- sd(sims_i)}

df <- df %>% mutate(Parameter = recode_factor(rownames(df),
                                              'sigma'='$\\sigma^{2}$',
                                              'sd_prov__Intercept'="$\\sigma^{2}_{P}$",
                                              'sd_prov:clon__Intercept'='$\\sigma^{2}_{G}$',
                                              'sd_prov:site__Intercept'='$\\sigma^{2}_{Inter}$',
                                              'sd_site__Intercept'='$\\sigma^{2}_{S}$',
                                              'sd_site:block__Intercept'='$\\sigma^{2}_{B}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)


pars <- mod %>% get_variables() %>%  as.vector() %>% str_subset("^(b)")
sims <- as.array(mod, pars = pars, fixed = TRUE)

df2 <- as.data.frame(matrix(NA,length(pars),4,dimnames = list(pars,c("Median","SD","InfCI","SupCI"))))
for(i in pars){
  sims_i <- sims[, , i]
  quan <- unname(quantile(sims_i, probs = probs))
  df2[i,"InfCI"] <- quan[1]
  df2[i,"SupCI"] <- quan[2]
  df2[i,"Median"] <- median(sims_i)
  df2[i,"SD"] <- sd(sims_i)}

df2 <- df2 %>% mutate(Parameter = recode_factor(rownames(df2),
                                                'b_age.sc'="$\\beta_{age}$",
                                                'b_Iage.scE2'= '$\\beta_{age2}$',
                                                'b_Intercept'='$\\beta_{0}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)


tab <- bind_rows(df,df2)

print(xtable(tab, type = "latex",digits=3), 
      file = paste0("tables/Posteriors/M2_MainVarPost.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})


# >> Figure S12. ####
# -------------  "


# only for the P1 partition


# >>>> Panel All sites ####

POST <- posterior_samples(mod,pars = "^r_prov\\[")
colnames(POST) <- str_sub(colnames(POST),8,-12)
POST <- as.data.frame(t(POST))
POST$prov <- as.factor(rownames(POST))

data <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv") %>%  dplyr::filter(P1=="train")
data <- droplevels(data)
ps <- data %>%
  group_by(prov) %>% 
  summarise_at(vars(paste0(rep("Q",6),1:6)), mean)
ps$max.Q.prov <- colnames(ps[,2:7])[apply(ps[,2:7],1,which.max)]
ps$prov <- as.factor(ps$prov)

posteriorsimpelmodellong <- inner_join(POST, ps[,c("prov","max.Q.prov")],by="prov") %>%  
  as_tibble() %>% 
  gather(key = "key", value = "value", -prov,-max.Q.prov) %>%
  group_by(prov) %>%
  dplyr::mutate(meanperprov = mean(value)) %>%
  ungroup()

pm2_all <- ggplot()+
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
  scale_y_discrete(labels=c("ALT"=parse(text = TeX("$P_{ALT}$")),
                            "ARM"=parse(text = TeX("$P_{ARM}$")),
                            "ARN"=parse(text = TeX("$P_{ARN}$")),
                            "BAY"=parse(text = TeX("$P_{BAY}$")),
                            "BON"=parse(text = TeX("$P_{BON}$")),
                            "CAD"=parse(text = TeX("$P_{CAD}$")),
                            "CAR"=parse(text = TeX("$P_{CAR}$")),
                            "CAS"=parse(text = TeX("$P_{CAS}$")),
                            "CEN"=parse(text = TeX("$P_{CEN}$")),
                            "COC"=parse(text = TeX("$P_{COC}$")),
                            "COM"=parse(text = TeX("$P_{COM}$")),
                            "CUE"=parse(text = TeX("$P_{CUE}$")),
                            "HOU"=parse(text = TeX("$P_{HOU}$")),
                            "LAM"=parse(text = TeX("$P_{LAM}$")),
                            "LEI"=parse(text = TeX("$P_{LEI}$")),
                            "MAD"=parse(text = TeX("$P_{MAD}$")),
                            "MIM"=parse(text = TeX("$P_{MIM}$")),
                            "OLB"=parse(text = TeX("$P_{OLB}$")),
                            "OLO"=parse(text = TeX("$P_{OLO}$")),
                            "ORI"=parse(text = TeX("$P_{ORI}$")),
                            "PET"=parse(text = TeX("$P_{PET}$")),
                            "PIA"=parse(text = TeX("$P_{PIA}$")),
                            "PIE"=parse(text = TeX("$P_{PIE}$")),
                            "PLE"=parse(text = TeX("$P_{PLE}$")),
                            "PUE"=parse(text = TeX("$P_{PUE}$")),
                            "QUA"=parse(text = TeX("$P_{QUA}$")),
                            "SAC"=parse(text = TeX("$P_{SAC}$")),
                            "SAL"=parse(text = TeX("$P_{SAL}$")),
                            "SEG"=parse(text = TeX("$P_{SEG}$")),
                            "SIE"=parse(text = TeX("$P_{SIE}$")),
                            "STJ"=parse(text = TeX("$P_{STJ}$")),
                            "TAM"=parse(text = TeX("$P_{TAM}$")),
                            "VAL"=parse(text = TeX("$P_{VAL}$")),
                            "VER"=parse(text = TeX("$P_{VER}$")))) +
  
  coord_cartesian(c(-0.35,0.3))+
  
  scale_discrete_manual("vline_color",
                        values = c("blue", "red", "blue", "black"), 
                        breaks = c(1, 2),
                        labels = c("2.5 and 97.5th percentiles", "Mean"),
                        name = NULL) +
  scale_fill_manual(values=c("orangered3",
                             "gold2","darkorchid3",
                             "navyblue",
                             "turquoise2",
                             "green3"), labels = c("Q1: Northern Africa", 
                                                   "Q2: Corsica",
                                                   "Q3: Central Spain",
                                                   "Q4: French Atlantic",
                                                   "Q5: Iberian Atlantic",
                                                   "Q6: South-eastern Spain"),name="Gene pools") +
  labs(title     = "All sites",
       y        = "", 
       x        = TeX("Intercepts P_{p}")
  ) +
  theme_bw() + 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14), 
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))



# >>>> Panel Portugal ####

POST <- posterior_samples(mod,pars = "^r_prov:site\\[.*portugal")
colnames(POST) <- str_sub(colnames(POST),13,-21)
POST <- as.data.frame(t(POST))
POST$prov <- as.factor(rownames(POST))

posteriorsimpelmodellong <- inner_join(POST, ps[,c("prov","max.Q.prov")],by="prov") %>%  
  as_tibble() %>% 
  gather(key = "key", value = "value", -prov,-max.Q.prov)%>%
  group_by(prov) %>%
  dplyr::mutate(meanperprov = mean(value))%>%
  ungroup()

pm2_portugal <- ggplot()+
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
  
  scale_y_discrete(labels=c("ALT"=parse(text = TeX("$P_{ALT}$")),
                            "ARM"=parse(text = TeX("$P_{ARM}$")),
                            "ARN"=parse(text = TeX("$P_{ARN}$")),
                            "BAY"=parse(text = TeX("$P_{BAY}$")),
                            "BON"=parse(text = TeX("$P_{BON}$")),
                            "CAD"=parse(text = TeX("$P_{CAD}$")),
                            "CAR"=parse(text = TeX("$P_{CAR}$")),
                            "CAS"=parse(text = TeX("$P_{CAS}$")),
                            "CEN"=parse(text = TeX("$P_{CEN}$")),
                            "COC"=parse(text = TeX("$P_{COC}$")),
                            "COM"=parse(text = TeX("$P_{COM}$")),
                            "CUE"=parse(text = TeX("$P_{CUE}$")),
                            "HOU"=parse(text = TeX("$P_{HOU}$")),
                            "LAM"=parse(text = TeX("$P_{LAM}$")),
                            "LEI"=parse(text = TeX("$P_{LEI}$")),
                            "MAD"=parse(text = TeX("$P_{MAD}$")),
                            "MIM"=parse(text = TeX("$P_{MIM}$")),
                            "OLB"=parse(text = TeX("$P_{OLB}$")),
                            "OLO"=parse(text = TeX("$P_{OLO}$")),
                            "ORI"=parse(text = TeX("$P_{ORI}$")),
                            "PET"=parse(text = TeX("$P_{PET}$")),
                            "PIA"=parse(text = TeX("$P_{PIA}$")),
                            "PIE"=parse(text = TeX("$P_{PIE}$")),
                            "PLE"=parse(text = TeX("$P_{PLE}$")),
                            "PUE"=parse(text = TeX("$P_{PUE}$")),
                            "QUA"=parse(text = TeX("$P_{QUA}$")),
                            "SAC"=parse(text = TeX("$P_{SAC}$")),
                            "SAL"=parse(text = TeX("$P_{SAL}$")),
                            "SEG"=parse(text = TeX("$P_{SEG}$")),
                            "SIE"=parse(text = TeX("$P_{SIE}$")),
                            "STJ"=parse(text = TeX("$P_{STJ}$")),
                            "TAM"=parse(text = TeX("$P_{TAM}$")),
                            "VAL"=parse(text = TeX("$P_{VAL}$")),
                            "VER"=parse(text = TeX("$P_{VER}$")))) +
  
  coord_cartesian(c(-0.35,0.3))+
  
  scale_discrete_manual("vline_color",
                        values = c("blue", "red", "blue", "black"), 
                        breaks = c(1, 2),
                        labels = c("2.5 and 97.5th percentiles", "Mean"),
                        name = NULL) +
  scale_fill_manual(values=c("orangered3",
                             "gold2","darkorchid3",
                             "navyblue",
                             "turquoise2",
                             "green3"), labels = c("Q1: Northern Africa", 
                                                   "Q2: Corsica",
                                                   "Q3: Central Spain",
                                                   "Q4: French Atlantic",
                                                   "Q5: Iberian Atlantic",
                                                   "Q6: South-eastern Spain"),name="Gene pools") +
  labs(title     = "Portugal",
       y        = "", 
       x        = TeX("Intercepts P_{p,Portugal}")
  ) +
  theme_bw() + 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14), 
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))



# >>>> Panel Caceres ####

POST <- posterior_samples(mod,pars = "^r_prov:site\\[.*caceres")
colnames(POST) <- str_sub(colnames(POST),13,-20)
POST <- as.data.frame(t(POST))
POST$prov <- as.factor(rownames(POST))

posteriorsimpelmodellong <- inner_join(POST, ps[,c("prov","max.Q.prov")],by="prov") %>%  
  as_tibble() %>% 
  gather(key = "key", value = "value", -prov,-max.Q.prov)%>%
  group_by(prov) %>%
  dplyr::mutate(meanperprov = mean(value))%>%
  ungroup()

pm2_caceres <- ggplot()+
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
  
  scale_y_discrete(labels=c("ALT"=parse(text = TeX("$P_{ALT}$")),
                            "ARM"=parse(text = TeX("$P_{ARM}$")),
                            "ARN"=parse(text = TeX("$P_{ARN}$")),
                            "BAY"=parse(text = TeX("$P_{BAY}$")),
                            "BON"=parse(text = TeX("$P_{BON}$")),
                            "CAD"=parse(text = TeX("$P_{CAD}$")),
                            "CAR"=parse(text = TeX("$P_{CAR}$")),
                            "CAS"=parse(text = TeX("$P_{CAS}$")),
                            "CEN"=parse(text = TeX("$P_{CEN}$")),
                            "COC"=parse(text = TeX("$P_{COC}$")),
                            "COM"=parse(text = TeX("$P_{COM}$")),
                            "CUE"=parse(text = TeX("$P_{CUE}$")),
                            "HOU"=parse(text = TeX("$P_{HOU}$")),
                            "LAM"=parse(text = TeX("$P_{LAM}$")),
                            "LEI"=parse(text = TeX("$P_{LEI}$")),
                            "MAD"=parse(text = TeX("$P_{MAD}$")),
                            "MIM"=parse(text = TeX("$P_{MIM}$")),
                            "OLB"=parse(text = TeX("$P_{OLB}$")),
                            "OLO"=parse(text = TeX("$P_{OLO}$")),
                            "ORI"=parse(text = TeX("$P_{ORI}$")),
                            "PET"=parse(text = TeX("$P_{PET}$")),
                            "PIA"=parse(text = TeX("$P_{PIA}$")),
                            "PIE"=parse(text = TeX("$P_{PIE}$")),
                            "PLE"=parse(text = TeX("$P_{PLE}$")),
                            "PUE"=parse(text = TeX("$P_{PUE}$")),
                            "QUA"=parse(text = TeX("$P_{QUA}$")),
                            "SAC"=parse(text = TeX("$P_{SAC}$")),
                            "SAL"=parse(text = TeX("$P_{SAL}$")),
                            "SEG"=parse(text = TeX("$P_{SEG}$")),
                            "SIE"=parse(text = TeX("$P_{SIE}$")),
                            "STJ"=parse(text = TeX("$P_{STJ}$")),
                            "TAM"=parse(text = TeX("$P_{TAM}$")),
                            "VAL"=parse(text = TeX("$P_{VAL}$")),
                            "VER"=parse(text = TeX("$P_{VER}$")))) +
  coord_cartesian(c(-0.35,0.3))+
  
  scale_discrete_manual("vline_color",
                        values = c("blue", "red", "blue", "black"), 
                        breaks = c(1, 2),
                        labels = c("2.5 and 97.5th percentiles", "Mean"),
                        name = NULL) +
  scale_fill_manual(values=c("orangered3",
                             "gold2","darkorchid3",
                             "navyblue",
                             "turquoise2",
                             "green3"), labels = c("Q1: Northern Africa", 
                                                   "Q2: Corsica",
                                                   "Q3: Central Spain",
                                                   "Q4: French Atlantic",
                                                   "Q5: Iberian Atlantic",
                                                   "Q6: South-eastern Spain"),name="Gene pools") +
  labs(title     = "Caceres",
       y        = "", 
       x        = TeX("Intercepts P_{p,Caceres}")
  ) +
  theme_bw() + 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14), 
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))



# >>>> Panel Madrid ####


POST <- posterior_samples(mod,pars = "^r_prov:site\\[.*madrid")
colnames(POST) <- str_sub(colnames(POST),13,-19)
POST <- as.data.frame(t(POST))
POST$prov <- as.factor(rownames(POST))

posteriorsimpelmodellong <- inner_join(POST, ps[,c("prov","max.Q.prov")],by="prov") %>%  
  as_tibble() %>% 
  gather(key = "key", value = "value", -prov,-max.Q.prov) %>%
  group_by(prov) %>%
  dplyr::mutate(meanperprov = mean(value)) %>%
  ungroup()

pm2_madrid <- ggplot()+
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
  
  coord_cartesian(c(-0.35,0.3))+
  
  scale_y_discrete(labels=c("ALT"=parse(text = TeX("$P_{ALT}$")),
                            "ARM"=parse(text = TeX("$P_{ARM}$")),
                            "ARN"=parse(text = TeX("$P_{ARN}$")),
                            "BAY"=parse(text = TeX("$P_{BAY}$")),
                            "BON"=parse(text = TeX("$P_{BON}$")),
                            "CAD"=parse(text = TeX("$P_{CAD}$")),
                            "CAR"=parse(text = TeX("$P_{CAR}$")),
                            "CAS"=parse(text = TeX("$P_{CAS}$")),
                            "CEN"=parse(text = TeX("$P_{CEN}$")),
                            "COC"=parse(text = TeX("$P_{COC}$")),
                            "COM"=parse(text = TeX("$P_{COM}$")),
                            "CUE"=parse(text = TeX("$P_{CUE}$")),
                            "HOU"=parse(text = TeX("$P_{HOU}$")),
                            "LAM"=parse(text = TeX("$P_{LAM}$")),
                            "LEI"=parse(text = TeX("$P_{LEI}$")),
                            "MAD"=parse(text = TeX("$P_{MAD}$")),
                            "MIM"=parse(text = TeX("$P_{MIM}$")),
                            "OLB"=parse(text = TeX("$P_{OLB}$")),
                            "OLO"=parse(text = TeX("$P_{OLO}$")),
                            "ORI"=parse(text = TeX("$P_{ORI}$")),
                            "PET"=parse(text = TeX("$P_{PET}$")),
                            "PIA"=parse(text = TeX("$P_{PIA}$")),
                            "PIE"=parse(text = TeX("$P_{PIE}$")),
                            "PLE"=parse(text = TeX("$P_{PLE}$")),
                            "PUE"=parse(text = TeX("$P_{PUE}$")),
                            "QUA"=parse(text = TeX("$P_{QUA}$")),
                            "SAC"=parse(text = TeX("$P_{SAC}$")),
                            "SAL"=parse(text = TeX("$P_{SAL}$")),
                            "SEG"=parse(text = TeX("$P_{SEG}$")),
                            "SIE"=parse(text = TeX("$P_{SIE}$")),
                            "STJ"=parse(text = TeX("$P_{STJ}$")),
                            "TAM"=parse(text = TeX("$P_{TAM}$")),
                            "VAL"=parse(text = TeX("$P_{VAL}$")),
                            "VER"=parse(text = TeX("$P_{VER}$")))) +
  scale_discrete_manual("vline_color",
                        values = c("blue", "red", "blue", "black"), 
                        breaks = c(1, 2),
                        labels = c("2.5 and 97.5th percentiles", "Mean"),
                        name = NULL) +
  scale_fill_manual(values=c("orangered3",
                             "gold2","darkorchid3",
                             "navyblue",
                             "turquoise2",
                             "green3"), labels = c("Q1: Northern Africa", 
                                                   "Q2: Corsica",
                                                   "Q3: Central Spain",
                                                   "Q4: French Atlantic",
                                                   "Q5: Iberian Atlantic",
                                                   "Q6: South-eastern Spain"),name="Gene pools") +
  labs(title     = "Madrid",
       y        = "", 
       x        = TeX("Intercepts P_{p,Madrid}")
  ) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14), 
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))



# >>>> Panel Bordeaux ####

POST <- posterior_samples(mod,pars = "^r_prov:site\\[.*bordeaux")
colnames(POST) <- str_sub(colnames(POST),13,-21)
POST <- as.data.frame(t(POST))
POST$prov <- as.factor(rownames(POST))

posteriorsimpelmodellong <- inner_join(POST, ps[,c("prov","max.Q.prov")],by="prov") %>%  
  as_tibble() %>% 
  gather(key = "key", value = "value", -prov,-max.Q.prov) %>%
  group_by(prov) %>%
  dplyr::mutate(meanperprov = mean(value)) %>%
  ungroup()

pm2_bordeaux <- ggplot() +
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
  
  coord_cartesian(c(-0.35,0.3))+
  
  scale_y_discrete(labels=c("ALT"=parse(text = TeX("$P_{ALT}$")),
                            "ARM"=parse(text = TeX("$P_{ARM}$")),
                            "ARN"=parse(text = TeX("$P_{ARN}$")),
                            "BAY"=parse(text = TeX("$P_{BAY}$")),
                            "BON"=parse(text = TeX("$P_{BON}$")),
                            "CAD"=parse(text = TeX("$P_{CAD}$")),
                            "CAR"=parse(text = TeX("$P_{CAR}$")),
                            "CAS"=parse(text = TeX("$P_{CAS}$")),
                            "CEN"=parse(text = TeX("$P_{CEN}$")),
                            "COC"=parse(text = TeX("$P_{COC}$")),
                            "COM"=parse(text = TeX("$P_{COM}$")),
                            "CUE"=parse(text = TeX("$P_{CUE}$")),
                            "HOU"=parse(text = TeX("$P_{HOU}$")),
                            "LAM"=parse(text = TeX("$P_{LAM}$")),
                            "LEI"=parse(text = TeX("$P_{LEI}$")),
                            "MAD"=parse(text = TeX("$P_{MAD}$")),
                            "MIM"=parse(text = TeX("$P_{MIM}$")),
                            "OLB"=parse(text = TeX("$P_{OLB}$")),
                            "OLO"=parse(text = TeX("$P_{OLO}$")),
                            "ORI"=parse(text = TeX("$P_{ORI}$")),
                            "PET"=parse(text = TeX("$P_{PET}$")),
                            "PIA"=parse(text = TeX("$P_{PIA}$")),
                            "PIE"=parse(text = TeX("$P_{PIE}$")),
                            "PLE"=parse(text = TeX("$P_{PLE}$")),
                            "PUE"=parse(text = TeX("$P_{PUE}$")),
                            "QUA"=parse(text = TeX("$P_{QUA}$")),
                            "SAC"=parse(text = TeX("$P_{SAC}$")),
                            "SAL"=parse(text = TeX("$P_{SAL}$")),
                            "SEG"=parse(text = TeX("$P_{SEG}$")),
                            "SIE"=parse(text = TeX("$P_{SIE}$")),
                            "STJ"=parse(text = TeX("$P_{STJ}$")),
                            "TAM"=parse(text = TeX("$P_{TAM}$")),
                            "VAL"=parse(text = TeX("$P_{VAL}$")),
                            "VER"=parse(text = TeX("$P_{VER}$")))) +
  scale_discrete_manual("vline_color",
                        values = c("blue", "red", "blue", "black"), 
                        breaks = c(1, 2),
                        labels = c("2.5 and 97.5th percentiles", "Mean"),
                        name = NULL) +
  scale_fill_manual(values=c("orangered3",
                             "gold2","darkorchid3",
                             "navyblue",
                             "turquoise2",
                             "green3"), labels = c("Q1: Northern Africa", 
                                                   "Q2: Corsica",
                                                   "Q3: Central Spain",
                                                   "Q4: French Atlantic",
                                                   "Q5: Iberian Atlantic",
                                                   "Q6: South-eastern Spain"),name="Gene pools") +
  labs(title = "Bordeaux",
       y        = "", 
       x        = TeX("Intercepts P_{p,Bordeaux}")
  ) +
  theme_bw() + 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14), 
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))


# >>>> Panel Asturias ####


POST <- posterior_samples(mod,pars = "^r_prov:site\\[.*asturias")
colnames(POST) <- str_sub(colnames(POST),13,-21)
POST <- as.data.frame(t(POST))
POST$prov <- as.factor(rownames(POST))

posteriorsimpelmodellong <- inner_join(POST, ps[,c("prov","max.Q.prov")],by="prov") %>%  
  as_tibble() %>% 
  gather(key = "key", value = "value", -prov,-max.Q.prov) %>%
  group_by(prov) %>%
  dplyr::mutate(meanperprov = mean(value)) %>%
  ungroup()

pm2_asturias <- ggplot() +
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
  
  coord_cartesian(c(-0.35,0.3))+
  
  scale_y_discrete(labels=c("ALT"=parse(text = TeX("$P_{ALT}$")),
                            "ARM"=parse(text = TeX("$P_{ARM}$")),
                            "ARN"=parse(text = TeX("$P_{ARN}$")),
                            "BAY"=parse(text = TeX("$P_{BAY}$")),
                            "BON"=parse(text = TeX("$P_{BON}$")),
                            "CAD"=parse(text = TeX("$P_{CAD}$")),
                            "CAR"=parse(text = TeX("$P_{CAR}$")),
                            "CAS"=parse(text = TeX("$P_{CAS}$")),
                            "CEN"=parse(text = TeX("$P_{CEN}$")),
                            "COC"=parse(text = TeX("$P_{COC}$")),
                            "COM"=parse(text = TeX("$P_{COM}$")),
                            "CUE"=parse(text = TeX("$P_{CUE}$")),
                            "HOU"=parse(text = TeX("$P_{HOU}$")),
                            "LAM"=parse(text = TeX("$P_{LAM}$")),
                            "LEI"=parse(text = TeX("$P_{LEI}$")),
                            "MAD"=parse(text = TeX("$P_{MAD}$")),
                            "MIM"=parse(text = TeX("$P_{MIM}$")),
                            "OLB"=parse(text = TeX("$P_{OLB}$")),
                            "OLO"=parse(text = TeX("$P_{OLO}$")),
                            "ORI"=parse(text = TeX("$P_{ORI}$")),
                            "PET"=parse(text = TeX("$P_{PET}$")),
                            "PIA"=parse(text = TeX("$P_{PIA}$")),
                            "PIE"=parse(text = TeX("$P_{PIE}$")),
                            "PLE"=parse(text = TeX("$P_{PLE}$")),
                            "PUE"=parse(text = TeX("$P_{PUE}$")),
                            "QUA"=parse(text = TeX("$P_{QUA}$")),
                            "SAC"=parse(text = TeX("$P_{SAC}$")),
                            "SAL"=parse(text = TeX("$P_{SAL}$")),
                            "SEG"=parse(text = TeX("$P_{SEG}$")),
                            "SIE"=parse(text = TeX("$P_{SIE}$")),
                            "STJ"=parse(text = TeX("$P_{STJ}$")),
                            "TAM"=parse(text = TeX("$P_{TAM}$")),
                            "VAL"=parse(text = TeX("$P_{VAL}$")),
                            "VER"=parse(text = TeX("$P_{VER}$")))) +
  scale_discrete_manual("vline_color",
                        values = c("blue", "red", "blue", "black"), 
                        breaks = c(1, 2),
                        labels = c("2.5 and 97.5th percentiles", "Mean"),
                        name = NULL) +
  scale_fill_manual(values=c("orangered3",
                             "gold2","darkorchid3",
                             "navyblue",
                             "turquoise2",
                             "green3"), labels = c("Q1: Northern Africa", 
                                                   "Q2: Corsica",
                                                   "Q3: Central Spain",
                                                   "Q4: French Atlantic",
                                                   "Q5: Iberian Atlantic",
                                                   "Q6: South-eastern Spain"),name="Gene pools") +
  labs(title = "Asturias",
       y        = "", 
       x        = TeX("Intercepts P_{p,Asturias}")
  ) +
  theme_bw() + 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14), 
        legend.text = element_text(size=18),
        legend.title = element_text(size=20))


# >>>> Merging the panels ####


pp2_all <- pm2_all + theme(legend.position = "none")
pp2_asturias <- pm2_asturias + theme(legend.position = "none")
pp2_bordeaux <- pm2_bordeaux + theme(legend.position = "none")
pp2_caceres <- pm2_caceres + theme(legend.position = "none")
pp2_portugal <- pm2_portugal + theme(legend.position = "none")
pp2_madrid <- pm2_madrid + theme(legend.position = "none")

g <-  ggarrange(pp2_all,pp2_asturias,pp2_bordeaux,pp2_caceres,pp2_portugal,pp2_madrid ,
                nrow=1)

ggsave(g,file="figs/SuppInfo/M2_SiteProvIntercepts.png",width=20,height=12) 


# Model M3 ####
# ======== "

# only for the P1 partition

mod <- readRDS(file="outputs/models/P1/MOD3.rds")


# >> Table S19. ####
# ------------- "

pars <- mod %>% get_variables() %>%  as.vector() %>% str_subset("^(sd|sigma)")
sims <- as.array(mod, pars = pars, fixed = TRUE)

df <- as.data.frame(matrix(NA,length(pars),4,dimnames = list(pars,c("Median","SD","InfCI","SupCI"))))
for(i in pars){
  sims_i <- square(sims[, , i])
  quan <- unname(quantile(sims_i, probs = probs))
  df[i,"InfCI"] <- quan[1]
  df[i,"SupCI"] <- quan[2]
  df[i,"Median"] <- median(sims_i)
  df[i,"SD"] <- sd(sims_i)}

df <- df %>% mutate(Parameter = recode_factor(rownames(df),
                                              'sigma'='$\\sigma^{2}$',
                                              'sd_prov__Intercept'="$\\sigma^{2}_{P}$",
                                              'sd_site_age__Intercept'='$\\sigma^{2}_{cs_{is}}$',
                                              'sd_clon__Intercept'='$\\sigma^{2}_{G}$',
                                              'sd_site__Intercept'='$\\sigma^{2}_{S}$',
                                              'sd_site:block__Intercept'='$\\sigma^{2}_{B}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)


pars <- mod %>% get_variables() %>%  as.vector() %>% str_subset("^(b)")
sims <- as.array(mod, pars = pars, fixed = TRUE)

df2 <- as.data.frame(matrix(NA,length(pars),4,dimnames = list(pars,c("Median","SD","InfCI","SupCI"))))
for(i in pars){
  sims_i <- sims[, , i]
  quan <- unname(quantile(sims_i, probs = probs))
  df2[i,"InfCI"] <- quan[1]
  df2[i,"SupCI"] <- quan[2]
  df2[i,"Median"] <- median(sims_i)
  df2[i,"SD"] <- sd(sims_i)}

df2 <- df2 %>% mutate(Parameter = recode_factor(rownames(df2),
                                                'b_age.sc'="$\\beta_{age}$",
                                                'b_Iage.scE2'= '$\\beta_{age2}$',
                                                'b_Intercept'='$\\beta_{0}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)

tab <- bind_rows(df,df2)

print(xtable(tab, type = "latex",digits=3), 
      file = paste0("tables/Posteriors/M3_MainVarPost.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})


# >> Table S20. ####
# ------------- "

df <- mod %>% 
  broom::tidyMCMC(estimate.method = "median",conf.int = T,conf.level = 0.95) %>% 
  filter(str_detect(term, "^(r_site\\[)")) %>% 
  rename_all(str_to_title) %>% 
  dplyr::rename("Median"=Estimate,"SD"=Std.error,"InfCI"=Conf.low,"SupCI"=Conf.high) %>% 
  mutate(Parameter = recode_factor(Term,
                                   'r_site[asturias,Intercept]'="$S_{Asturias}$",
                                   'r_site[bordeaux,Intercept]'= '$S_{Bordeaux}$',
                                   'r_site[caceres,Intercept]'='$S_{Caceres}$',
                                   'r_site[madrid,Intercept]'='$S_{Madrid}$',
                                   'r_site[portugal,Intercept]'="$S_{Portugal}$")) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)

print(xtable(df, type = "latex",digits=3), 
      file = paste0("tables/Posteriors/M3_SiteInterceptsPost.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})


# >> Table S21. ####
# ------------- "


df <- mod %>% 
  broom::tidyMCMC(estimate.method = "median",conf.int = T,conf.level = 0.95) %>% 
  filter(str_detect(term, "^(r_site_age\\[)")) %>% 
  rename_all(str_to_title) %>% 
  dplyr::rename("Median"=Estimate,"SD"=Std.error,"InfCI"=Conf.low,"SupCI"=Conf.high) %>% 
  mutate(Parameter = recode_factor(Term,
                                   'r_site_age[asturias10,Intercept]' = '$cs_{1,Asturias}$',
                                   'r_site_age[portugal27,Intercept]' = '$cs_{4,Portugal}$',
                                   'r_site_age[portugal20,Intercept]' = '$cs_{3,Portugal}$',
                                   'r_site_age[asturias21,Intercept]' = '$cs_{2,Asturias}$',
                                   'r_site_age[portugal11,Intercept]' = '$cs_{1,Portugal}$',
                                   'r_site_age[madrid13,Intercept]' = '$cs_{1,Madrid}$',
                                   'r_site_age[asturias37,Intercept]' = '$cs_{3,Asturias}$',
                                   'r_site_age[portugal15,Intercept]' = '$cs_{2,Portugal}$',
                                   'r_site_age[bordeaux25,Intercept]' = '$cs_{1,Bordeaux}$',
                                   'r_site_age[bordeaux37,Intercept]' = '$cs_{2,Bordeaux}$',
                                   'r_site_age[caceres8,Intercept]' = '$cs_{1,Caceres}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)

print(xtable(df, type = "latex",digits=3), 
      file = paste0("tables/Posteriors/M3_SiteClimSimInterceptsPost.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})



# Model M3bis ####
# =========== "

# only for the P1 partition

mod <- readRDS(file="outputs/models/P1/MOD13.rds")


# >> Table S23. ####
# ------------- "

pars <- mod %>% get_variables() %>%  as.vector() %>% str_subset("^(sd|sigma)")
sims <- as.array(mod, pars = pars, fixed = TRUE)

df <- as.data.frame(matrix(NA,length(pars),4,dimnames = list(pars,c("Median","SD","InfCI","SupCI"))))
for(i in pars){
  sims_i <- square(sims[, , i])
  quan <- unname(quantile(sims_i, probs = probs))
  df[i,"InfCI"] <- quan[1]
  df[i,"SupCI"] <- quan[2]
  df[i,"Median"] <- median(sims_i)
  df[i,"SD"] <- sd(sims_i)}

df <- df %>% mutate(Parameter = recode_factor(rownames(df),
                                              'sigma'='$\\sigma^{2}$',
                                              'sd_prov__Intercept'="$\\sigma^{2}_{P}$",
                                              'sd_site_age__Intercept'='$\\sigma^{2}_{cs_{is}}$',
                                              'sd_prov:clon__Intercept'='$\\sigma^{2}_{G}$',
                                              'sd_block__Intercept'='$\\sigma^{2}_{B}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)


pars <- mod %>% get_variables() %>%  as.vector() %>% str_subset("^(b)")
sims <- as.array(mod, pars = pars, fixed = TRUE)

df2 <- as.data.frame(matrix(NA,length(pars),4,dimnames = list(pars,c("Median","SD","InfCI","SupCI"))))

for(i in pars){
  sims_i <- sims[, , i]
  quan <- unname(quantile(sims_i, probs = probs))
  df2[i,"InfCI"] <- quan[1]
  df2[i,"SupCI"] <- quan[2]
  df2[i,"Median"] <- median(sims_i)
  df2[i,"SD"] <- sd(sims_i)}

df2 <- df2 %>% mutate(Parameter = recode_factor(rownames(df2),
                                                'b_age.sc'="$\\beta_{age}$",
                                                'b_Iage.scE2'= '$\\beta_{age2}$',
                                                'b_Intercept'='$\\beta_{0}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)

tab <- bind_rows(df,df2)

print(xtable(tab, type = "latex",digits=3), 
      file = paste0("tables/Posteriors/M3bis_MainVarPost.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})


# >> Table S24. ####
# ------------- "

df <- mod %>% 
  broom::tidyMCMC(estimate.method = "median",conf.int = T,conf.level = 0.95) %>% 
  filter(str_detect(term, "^(r_site_age\\[)")) %>% 
  rename_all(str_to_title) %>% 
  dplyr::rename("Median"=Estimate,"SD"=Std.error,"InfCI"=Conf.low,"SupCI"=Conf.high) %>% 
  mutate(Parameter = recode_factor(Term,
                                   'r_site_age[asturias10,Intercept]' = '$cs_{1,Asturias}$',
                                   'r_site_age[portugal27,Intercept]' = '$cs_{4,Portugal}$',
                                   'r_site_age[portugal20,Intercept]' = '$cs_{3,Portugal}$',
                                   'r_site_age[asturias21,Intercept]' = '$cs_{2,Asturias}$',
                                   'r_site_age[portugal11,Intercept]' = '$cs_{1,Portugal}$',
                                   'r_site_age[madrid13,Intercept]' = '$cs_{1,Madrid}$',
                                   'r_site_age[asturias37,Intercept]' = '$cs_{3,Asturias}$',
                                   'r_site_age[portugal15,Intercept]' = '$cs_{2,Portugal}$',
                                   'r_site_age[bordeaux25,Intercept]' = '$cs_{1,Bordeaux}$',
                                   'r_site_age[bordeaux37,Intercept]' = '$cs_{2,Bordeaux}$',
                                   'r_site_age[caceres8,Intercept]' = '$cs_{1,Caceres}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)

print(xtable(df, type = "latex",digits=3), 
      file = paste0("tables/Posteriors/M3bis_SiteClimSimInterceptsPost.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})



# Model M4 ####
# ======== "

# only for the P1 partition

mod <- readRDS(file="outputs/models/P1/MOD4.rds")


# >> Table S25. ####
# ------------- "

pars <- mod %>% get_variables() %>%  as.vector() %>% str_subset("^(sd|sigma)")
sims <- as.array(mod, pars = pars, fixed = TRUE)

df <- as.data.frame(matrix(NA,length(pars),4,dimnames = list(pars,c("Median","SD","InfCI","SupCI"))))
for(i in pars){
  sims_i <- square(sims[, , i])
  quan <- unname(quantile(sims_i, probs = probs))
  df[i,"InfCI"] <- quan[1]
  df[i,"SupCI"] <- quan[2]
  df[i,"Median"] <- median(sims_i)
  df[i,"SD"] <- sd(sims_i)}

df <- df %>% mutate(Parameter = recode_factor(rownames(df),
                                              'sigma'='$\\sigma^{2}$',
                                              'sd_prov__Intercept'="$\\sigma^{2}_{P}$",
                                              'sd_site_age__Intercept'='$\\sigma^{2}_{cs_{is}}$',
                                              'sd_mmQ1Q2Q3Q4Q5Q6__Intercept'='$\\sigma^{2}_{g_{j}}$',
                                              'sd_prov:clon__Intercept'='$\\sigma^{2}_{G}$',
                                              'sd_site__Intercept'='$\\sigma^{2}_{S}$',
                                              'sd_site:block__Intercept'='$\\sigma^{2}_{B}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)


pars <- mod %>% get_variables() %>%  as.vector() %>% str_subset("^(b)")
sims <- as.array(mod, pars = pars, fixed = TRUE)

df2 <- as.data.frame(matrix(NA,length(pars),4,dimnames = list(pars,c("Median","SD","InfCI","SupCI"))))

for(i in pars){
  sims_i <- sims[, , i]
  quan <- unname(quantile(sims_i, probs = probs))
  df2[i,"InfCI"] <- quan[1]
  df2[i,"SupCI"] <- quan[2]
  df2[i,"Median"] <- median(sims_i)
  df2[i,"SD"] <- sd(sims_i)}

df2 <- df2 %>% mutate(Parameter = recode_factor(rownames(df2),
                                                'b_age.sc'="$\\beta_{age}$",
                                                'b_Iage.scE2'= '$\\beta_{age2}$',
                                                'b_Intercept'='$\\beta_{0}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)


tab <- bind_rows(df,df2)

print(xtable(tab, type = "latex",digits=3), 
      file = paste0("tables/Posteriors/M4_MainVarPost.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})



# Model M5 ####
# ======== "

# only for the P1 partition

mod <- readRDS(file="outputs/models/P1/MOD5.rds")


# >> Table S26. ####
# ------------- "

pars <- mod %>% get_variables() %>%  as.vector() %>% str_subset("^(sd|sigma)")
sims <- as.array(mod, pars = pars, fixed = TRUE)

df <- as.data.frame(matrix(NA,length(pars),4,dimnames = list(pars,c("Median","SD","InfCI","SupCI"))))

for(i in pars){
  sims_i <- square(sims[, , i])
  quan <- unname(quantile(sims_i, probs = probs))
  df[i,"InfCI"] <- quan[1]
  df[i,"SupCI"] <- quan[2]
  df[i,"Median"] <- median(sims_i)
  df[i,"SD"] <- sd(sims_i)}

df <- df %>% mutate(Parameter = recode_factor(rownames(df),
                                              'sigma'='$\\sigma^{2}$',
                                              'sd_prov__Intercept'="$\\sigma^{2}_{P}$",
                                              'sd_site_age__Intercept'='$\\sigma^{2}_{cs_{is}}$',
                                              'sd_clon1__Intercept'='$\\sigma^{2}_{A_{NA}}$',
                                              'sd_clon2__Intercept'='$\\sigma^{2}_{A_{C}}$',
                                              'sd_clon3__Intercept'='$\\sigma^{2}_{A_{CS}}$',
                                              'sd_clon4__Intercept'='$\\sigma^{2}_{A_{FA}}$',
                                              'sd_clon5__Intercept'='$\\sigma^{2}_{A_{IA}}$',
                                              'sd_clon6__Intercept'='$\\sigma^{2}_{A_{SES}}$',
                                              'sd_mmQ1Q2Q3Q4Q5Q6__Intercept'='$\\sigma^{2}_{g_{j}}$',
                                              'sd_prov:clon__Intercept'='$\\sigma^{2}_{G}$',
                                              'sd_site__Intercept'='$\\sigma^{2}_{S}$',
                                              'sd_site:block__Intercept'='$\\sigma^{2}_{B}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)


pars <- mod %>% get_variables() %>%  as.vector() %>% str_subset("^(b)")
sims <- as.array(mod, pars = pars, fixed = TRUE)

df2 <- as.data.frame(matrix(NA,length(pars),4,dimnames = list(pars,c("Median","SD","InfCI","SupCI"))))

for(i in pars){
  sims_i <- sims[, , i]
  quan <- unname(quantile(sims_i, probs = probs))
  df2[i,"InfCI"] <- quan[1]
  df2[i,"SupCI"] <- quan[2]
  df2[i,"Median"] <- median(sims_i)
  df2[i,"SD"] <- sd(sims_i)}

df2 <- df2 %>% mutate(Parameter = recode_factor(rownames(df2),
                                                'b_age.sc'="$\\beta_{age}$",
                                                'b_Iage.scE2'= '$\\beta_{age2}$',
                                                'b_Intercept'='$\\beta_{0}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)

tab <- bind_rows(df,df2)

print(xtable(tab, type = "latex",digits=3), 
      file = paste0("tables/Posteriors/M5_MainVarPost.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})

# Model M6 ####
# ======== "

# only for the P1 partition

mod <- readRDS(file="outputs/models/P1/MOD6.rds")


# >> Table S29. ####
# ------------- "

pars <- mod %>% get_variables() %>%  as.vector() %>% str_subset("^(sd|sigma)")
sims <- as.array(mod, pars = pars, fixed = TRUE)

df <- as.data.frame(matrix(NA,length(pars),4,dimnames = list(pars,c("Median","SD","InfCI","SupCI"))))

for(i in pars){
  sims_i <- square(sims[, , i])
  quan <- unname(quantile(sims_i, probs = probs))
  df[i,"InfCI"] <- quan[1]
  df[i,"SupCI"] <- quan[2]
  df[i,"Median"] <- median(sims_i)
  df[i,"SD"] <- sd(sims_i)}

df <- df %>% mutate(Parameter = recode_factor(rownames(df),
                                              'sigma'='$\\sigma^{2}$',
                                              'sd_prov__Intercept'="$\\sigma^{2}_{P}$",
                                              'sd_site_age__Intercept'='$\\sigma^{2}_{cs_{is}}$',
                                              'sd_prov_clim__Intercept'='$\\sigma^{2}_{cp_{p}}$',
                                              'sd_mmQ1Q2Q3Q4Q5Q6__Intercept'='$\\sigma^{2}_{g_{j}}$',
                                              'sd_prov:clon__Intercept'='$\\sigma^{2}_{G}$',
                                              'sd_site__Intercept'='$\\sigma^{2}_{S}$',
                                              'sd_site:block__Intercept'='$\\sigma^{2}_{B}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)


pars <- mod %>% get_variables() %>%  as.vector() %>% str_subset("^(b)")
sims <- as.array(mod, pars = pars, fixed = TRUE)

df2 <- as.data.frame(matrix(NA,length(pars),4,dimnames = list(pars,c("Median","SD","InfCI","SupCI"))))

for(i in pars){
  sims_i <- sims[, , i]
  quan <- unname(quantile(sims_i, probs = probs))
  df2[i,"InfCI"] <- quan[1]
  df2[i,"SupCI"] <- quan[2]
  df2[i,"Median"] <- median(sims_i)
  df2[i,"SD"] <- sd(sims_i)}

df2 <- df2 %>% mutate(Parameter = recode_factor(rownames(df2),
                                                'b_age.sc'="$\\beta_{age}$",
                                                'b_Iage.scE2'= '$\\beta_{age2}$',
                                                'b_Intercept'='$\\beta_{0}$')) %>% 
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)

tab <- bind_rows(df,df2)

print(xtable(tab, type = "latex",digits=3), 
      file = paste0("tables/Posteriors/M6_MainVarPost.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})

# Model M7 ####
# ======== "
# Model M8 ####
# ======== "
# Model M9 ####
# ======== "
# Model M10 ####
# ========= "
# Model M11 ####
# ========= "
# Model M12 ####
# ========= "


