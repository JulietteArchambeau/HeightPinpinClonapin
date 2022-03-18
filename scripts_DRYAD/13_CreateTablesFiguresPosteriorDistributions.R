########################################################################################################################"
#                                                                                                                      #
#               Tables and figures of the posterior distributions in the Supplementary Information                     #
#                                                                                                                      #
#                                   Juliette Archambeau                                                                #
#                                       18/03/2022                                                                     #
#                                                                                                                      #
########################################################################################################################"

library(broom)     # CRAN v0.5.2
library(latex2exp) # CRAN v0.4.0
library(ggplot2)   # CRAN v3.3.1
library(ggpubr)    # CRAN v0.2.1
library(tidybayes) # CRAN v2.0.1
library(dplyr)     # CRAN v1.0.0
library(bayesplot) # CRAN v1.7.1
color_scheme_set("green")
library(xtable)    # CRAN v1.8-4
library(ggridges)  # CRAN v0.5.1
library(tidyverse) # CRAN v1.3.0
library(tibble)    # CRAN v2.1.3
library(brms)      # CRAN v2.11.1


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

# only for the P1 and P2 partition.

mod <- readRDS(file= paste0("outputs/models/",part,"/MOD7.rds"))

# >> Tables S30 and S39. ####
# --------------------- "

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
                                              'sd_mmQ1Q2Q3Q4Q5Q6__Intercept'='$\\sigma^{2}_{g_{j}}$',
                                              'sd_site__bio5_prov.sc'="$\\sigma^{2}_{\\beta_{max.temp,s}}$",
                                              'sd_site__bio14_prov.sc'="$\\sigma^{2}_{\\beta_{min.pre,s}}$",
                                              'sd_site__gPEA.sc'="$\\sigma^{2}_{\\beta_{gPEA,s}}$",
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
      file = paste0("tables/Posteriors/M7_MainVarPost.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})


# Model M8 ####
# ======== "


# only for the P1 and P2 partition.

mod <- readRDS(file= paste0("outputs/models/",part,"/MOD8.rds"))

# >> Tables S31 and S40. ####
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
                                              'sd_mmQ1Q2Q3Q4Q5Q6__Intercept'='$\\sigma^{2}_{g_{j}}$',
                                              'sd_site__bio5_prov.sc'="$\\sigma^{2}_{\\beta_{max.temp,s}}$",
                                              'sd_site__bio14_prov.sc'="$\\sigma^{2}_{\\beta_{min.pre,s}}$",
                                              'sd_site__rPEA.sc'="$\\sigma^{2}_{\\beta_{rPEA,s}}$",
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
      file = paste0("tables/Posteriors/M8_MainVarPost.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})



# >> Figures S15, S17 and S19 ####
# --------------------------- "

# for the three partitions

models <- list()
models[[1]] <- readRDS(file=paste0("outputs/models/",part,"/MOD7.rds"))
models[[2]] <- readRDS(file=paste0("outputs/models/",part,"/MOD8.rds"))
names(models) <- c("M7","M8")


figs.GP <- list()
figs.beta <- list()


# Panels a. Population structure 

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
         x = "") + 
    theme_bw() + theme(axis.text = element_text(size=22),
                       axis.title = element_text(size=22), 
                       legend.text = element_text(size=16),
                       legend.title = element_text(size=17))
}

figs.GP[[1]] <- figs.GP[[1]] + theme(legend.position = "none")
pGP <- ggarrange(figs.GP[[1]],figs.GP[[2]],
                 labels=c("M7. a)","M8. a)"),
                 font.label = list(size = 20),
                 hjust=c(-0.1,-0.1),
                 vjust=c(1.6,1.6),
                 nrow=1,
                 widths = c(1,1.3))


# Panels b. Provenance climates and PEAs

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
         x = "") + 
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
    
    theme_bw() + theme(axis.text = element_text(size=22),
                       axis.title = element_text(size=22),
                       legend.position = "none",
                       legend.text = element_text(size=16),
                       legend.title = element_text(size=17))
}


pbeta <- ggarrange(figs.beta[[1]],figs.beta[[2]],
                   labels=c("M7. b)","M8. b)"),
                   font.label = list(size = 20),
                   hjust=c(-0.1,-0.1),
                   vjust=c(1.6,1.6),
                   nrow=1)

# Merge the panels:
figtot <- ggarrange(pGP,pbeta,nrow=2,heights=c(1,2)) 

# Save the figure:
ggsave(figtot,file=paste0("figs/SuppInfo/M7M8Posteriors",part,".png"),height=12,width=20)




# Model M9 ####
# ======== "

# only for the P1 and P2 partition.

mod <- readRDS(file= paste0("outputs/models/",part,"/MOD9.rds"))


# >> Tables S32 and S41. ####
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
                                              'sd_mmQ1Q2Q3Q4Q5Q6__Intercept'='$\\sigma^{2}_{g_{j}}$',
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
      file = paste0("tables/Posteriors/M9_MainVarPost.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})


# Model M10 ####
# ========= "

# only for the P1 and P2 partition.

mod <- readRDS(file= paste0("outputs/models/",part,"/MOD10.rds"))


# >> Tables S33 and S42. ####
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
                                              'sd_site__Intercept'='$\\sigma^{2}_{S}$',
                                              'sd_site__bio5_prov.sc'="$\\sigma^{2}_{\\beta_{max.temp,s}}$",
                                              'sd_site__bio14_prov.sc'="$\\sigma^{2}_{\\beta_{min.pre,s}}$",
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
      file = paste0("tables/Posteriors/M10_MainVarPost.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})

# Model M11 ####
# ========= "

# only for the P1 and P2 partition.

mod <- readRDS(file= paste0("outputs/models/",part,"/MOD11.rds"))


# >> Tables S34 and S43. ####
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
                                              'sd_site__Intercept'='$\\sigma^{2}_{S}$',
                                              'sd_site__gPEA.sc'="$\\sigma^{2}_{\\beta_{gPEA,s}}$",
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
      file = paste0("tables/Posteriors/M11_MainVarPost.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})


# Model M12 ####
# ========= "

# only for the P1 and P2 partition.

mod <- readRDS(file= paste0("outputs/models/",part,"/MOD12.rds"))


# >> Tables S35 and S44. ####
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
                                              'sd_site__Intercept'='$\\sigma^{2}_{S}$',
                                              'sd_site__rPEA.sc'="$\\sigma^{2}_{\\beta_{rPEA,s}}$",
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
      file = paste0("tables/Posteriors/M12_MainVarPost.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})



# >> Figures S16, S18 and S20 ####
# --------------------------- "

# for the three partitions


# >>>> Panel M9 => Population structure

mod <- readRDS(file=paste0("outputs/models/",part,"/MOD9.rds"))

POST <- posterior_samples(mod,pars = "^r_mmQ1Q2Q3Q4Q5Q6\\[")
colnames(POST) <- str_sub(colnames(POST),18,-12)
POST <- as.data.frame(t(POST))
POST$genepool <- as.factor(rownames(POST))

posteriorsimpelmodellong <- POST %>%  
  as_tibble() %>% 
  gather(key = "key", value = "value", -genepool)%>%
  group_by(genepool) %>%
  dplyr::mutate(meanpergenepool = mean(value)) %>%
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
       x="") + 
  theme_bw() + theme(axis.text = element_text(size=18),
                     axis.title = element_text(size=18), 
                     legend.text = element_text(size=16),
                     legend.title = element_text(size=17))


# >>>> Panel M10 => Climate in the provenances

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

posteriorsimpelmodellong <- POST %>%  
  as_tibble() %>% 
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
  
  theme_bw() + theme(axis.text = element_text(size=22),
                     axis.title = element_text(size=22), 
                     legend.position = "none", 
                     plot.title = element_text(size=22)) 


# >>>> Panel M11 => gPEAs

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

posteriorsimpelmodellong <- POST %>%  
  as_tibble() %>% 
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
  theme_bw() + theme(axis.text = element_text(size=22),
                     axis.title = element_text(size=22), 
                     legend.text = element_text(size=16),
                     legend.title = element_text(size=17),
                     plot.title = element_text(size=22))  + 
  guides(vline_color = FALSE)


# >>>> Panel M12 => rPEAs


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

posteriorsimpelmodellong <- POST %>%  
  as_tibble() %>% 
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
       x="") + 
  scale_fill_manual(values=c(vir_lite("cyan2",ds=ds),
                             vir_lite("navyblue",ds=ds),
                             vir_lite("pink",ds=ds),
                             vir_lite("deeppink",ds=ds),
                             vir_lite("dodgerblue2",ds=ds))) +
  theme_bw() + theme(axis.text = element_text(size=22),
                     axis.title = element_text(size=22), 
                     legend.position = "none", 
                     plot.title = element_text(size=22)) 


p1 <- ggarrange(pGP,pCP,labels=c("M9","M10"),font.label = list(size = 20),nrow=1,widths = c(1.2,1))
p2 <- ggarrange(pgpea,prpea,labels=c("M11","M12"),font.label = list(size = 20),nrow=1,widths = c(1.2,1))
fig <- ggarrange(p1,p2,nrow=2)

if(part=="P1"){
  ggsave(fig, file=paste0("figs/manuscript/",part,"M9toM12PosteriorDistri.png"),width=20,height=12)  
} else{
  ggsave(fig, file=paste0("figs/SuppInfo/",part,"M9toM12PosteriorDistri.png"),width=20,height=12)  
}
