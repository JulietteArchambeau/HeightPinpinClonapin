##############################################################################################################"
##################                                                                    ########################"
##################     Fig. to understand gen and env components after review         ########################"
##################                                                                    ########################"
##############################################################################################################"

# This figure aims at showing the variance partitioning and understanding the drivers behind the genetic and plastic components.

library(dplyr)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(stringr)
library(rstanarm)
library(ggforce)
library(tidybayes)
library(latex2exp)
library(ggplot2)
library(brms)


## Variance partitioning

name.models <- c("MOD0","MOD1","MOD2")
models <- lapply(name.models,function(x) readRDS(file=paste0("outputs/models/P1/",x,".rds")))
names(models)  <- name.models

brms::bayes_R2(m2,re.form=NA)
brms::bayes_R2(m2,re_formula=~(1|site))
brms::bayes_R2(m2,re_formula=~(1|site/block))
brms::bayes_R2(m2,re_formula=~(1|prov/clon))
brms::bayes_R2(m2,re_formula=~(1|prov:site))


### Function to extract R2 without age
VarPart <- function(fit) {
  
  # Phenotypic variance
  y <- rstanarm::get_y(fit)
  var_y <- var(y)
  
  # Explained variance over the phenotypic variance
  ypred <- pp_expect(fit, transform = TRUE)
  var_ypred <- apply(ypred, 1, var)
  tab <- posterior_summary(var_ypred/var_y,probs = c(0.025, 0.975)) %>% 
    as_tibble() %>% 
    dplyr::mutate(param="R2all")
  
  # Variance explained by age over the phenotypic variance
  ypred <- pp_expect(fit, transform = TRUE,re.form=NA)
  var_ypred_age <- apply(ypred, 1, var)
  tab <- posterior_summary(var_ypred_age/var_y,probs = c(0.025, 0.975)) %>% 
    as_tibble() %>% 
    dplyr::mutate(param="R2fix") %>% 
    bind_rows(tab)
  
  # Variance explained by site/block over the phenotypic variance - variance explained by age
  ypred <- pp_expect(fit, transform = TRUE,re_formula=~(1|site/block))
  var_ypred <- apply(ypred, 1, var)
  tab <- posterior_summary((var_ypred-var_ypred_age)/(var_y-var_ypred_age),probs = c(0.025, 0.975)) %>% 
    as_tibble() %>% 
    dplyr::mutate(param="varEnv") %>% 
    bind_rows(tab)
  
  # Variance explained by prov/clon over the phenotypic variance - variance explained by age
  ypred <- pp_expect(fit, transform = TRUE,re_formula=~(1|prov/clon))
  var_ypred <- apply(ypred, 1, var)
  tab <- posterior_summary((var_ypred-var_ypred_age)/(var_y-var_ypred_age),probs = c(0.025, 0.975)) %>% 
    as_tibble() %>% 
    dplyr::mutate(param="varGen") %>% 
    bind_rows(tab)
  
  # Variance explained by prov:site over the phenotypic variance - variance explained by age
  ypred <- pp_expect(fit, transform = TRUE,re_formula=~(1|prov:site))
  var_ypred <- apply(ypred, 1, var)
  tab <- posterior_summary((var_ypred-var_ypred_age)/(var_y-var_ypred_age),probs = c(0.025, 0.975)) %>% 
    as_tibble() %>% 
    dplyr::mutate(param="varGxE") %>% 
    bind_rows(tab)
  
  # Residual and total explained variance over the phenotypic variance - variance explained by age 
  ypred <- pp_expect(fit, transform = TRUE)
  e <- -1 * sweep(ypred, 2, y)
  var_ypred <- apply(ypred, 1, var)
  var_e <- apply(e, 1, var)
  tab <- posterior_summary((var_ypred-var_ypred_age)/(var_y-var_ypred_age),probs = c(0.025, 0.975)) %>% 
    as_tibble() %>% 
    dplyr::mutate(param="varTot") %>% 
    bind_rows(tab)
  tab <- posterior_summary((var_e)/(var_y-var_ypred_age),probs = c(0.025, 0.975)) %>% 
    as_tibble() %>% 
    dplyr::mutate(param="varRes") %>% 
    bind_rows(tab)
}


tab <- lapply(models, VarPart)

tab <- VarPart(m1)

tabsub <- tab$MOD1 %>% dplyr::filter(param %in% c("varRes","varEnv","varGen")) 
tabsub$param <- as.factor(tabsub$param)
tabsub$param <- factor(tabsub$param,levels=c("varRes","varEnv","varGen"))

VarPart <- ggplot(tabsub) +
  geom_col(aes(x = param, y = Estimate, fill = param),alpha=0.8) + 
  scale_fill_manual(values=c("#EF8A62", "#7FBC41","#67A9CF")) + 
  xlab("") + ylab("Proportion of variance") +
  scale_x_discrete(labels=c("Residual","Environment","Genetic")) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size=25),
        axis.text = element_text(size = 25))


##### Environmental component

m3 <- readRDS(file="outputs/models/P1/MOD3.rds")

df <- m3 %>% broom::tidyMCMC(estimate.method = "median",conf.int = T,conf.level = 0.95) %>% 
  filter(str_detect(term, "^(r_site\\[)|^(r_site_age\\[)")) %>% 
  mutate(Param = case_when(
      startsWith(term, "r_site[") ~ "Intercepts of the common gardens",
      startsWith(term, "r_site_age") ~ "Intercepts of the annual climate in the common gardens"),
    Site=case_when(
      grepl("bordeaux",term) ~ "Bordeaux",
      grepl("madrid",term) ~ "Madrid",
      grepl("caceres",term)~ "Caceres",
      grepl("portugal",term)~ "Portugal",
      grepl("asturias",term)~ "Asturias"))

df$Param <- as.factor(df$Param)
df$Param <- factor(df$Param,levels=c("Intercepts of the common gardens","Intercepts of the annual climate in the common gardens"))

# Parameters of the funtion vir_lite
source("scripts/Functions/vir_lite.R")
ds=0.7


p <- df %>% ggplot(aes(x = term, y = estimate,ymin = conf.low, ymax = conf.high,color=Site)) +
  coord_flip() +
  facet_grid(.~Param,scales="free_y", space = "free") +
  #facet_grid(rows = vars(Param),scales="free_y", space = "free") +
  #facet_wrap(.~Param,scales="free_y",ncol=1) +
  geom_pointinterval(position = position_dodge(width = .6),point_size=2.5,size=2.5) +
  ylab("") + xlab("") +
  scale_color_manual(values=c(vir_lite("dodgerblue2",ds=ds),
                             vir_lite("deeppink",ds=ds),
                             vir_lite("pink",ds=ds),
                             vir_lite("navyblue",ds=ds),
                             vir_lite("cyan2",ds=ds))) +
  theme_bw() +
  scale_x_discrete(labels=c("r_site_age[caceres8,Intercept]"=parse(text = TeX("$cs_{8 months,Caceres}$")),
                            "r_site_age[bordeaux25,Intercept]"=parse(text = TeX("$cs_{25months,Bordeaux}$")),
                            "r_site_age[bordeaux37,Intercept]"=parse(text = TeX("$cs_{37months,Bordeaux}$")),
                            "r_site_age[portugal11,Intercept]"=parse(text = TeX("$cs_{11months,Portugal}$")),
                            "r_site_age[portugal15,Intercept]"=parse(text = TeX("$cs_{15months,Portugal}$")),
                            "r_site_age[portugal20,Intercept]"=parse(text = TeX("$cs_{20months,Portugal}$")),
                            "r_site_age[portugal27,Intercept]"=parse(text = TeX("$cs_{27months,Portugal}$")),
                            "r_site_age[madrid13,Intercept]"=parse(text = TeX("$cs_{13months,Madrid}$")),
                            "r_site_age[asturias10,Intercept]"=parse(text = TeX("$cs_{10months,Asturias}$")),
                            "r_site_age[asturias21,Intercept]"=parse(text = TeX("$cs_{21months,Asturias}$")),
                            "r_site_age[asturias37,Intercept]"=parse(text = TeX("$cs_{37months,Asturias}$")),
                            'r_site[asturias,Intercept]'=parse(text = TeX("$S_{Asturias}$")),
                            'r_site[bordeaux,Intercept]'=parse(text = TeX("$S_{Bordeaux}$")),
                            'r_site[caceres,Intercept]'=parse(text = TeX("$S_{Caceres}$")),
                            'r_site[madrid,Intercept]'=parse(text = TeX("$S_{Madrid}$")),
                            'r_site[portugal,Intercept]'=parse(text = TeX("$S_{Portugal}$")))) +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=16),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=11),
        legend.position = c(0.88,0.50),
        legend.background = element_rect(colour = "grey"),
        panel.grid.minor.y=element_blank(),
        strip.text.x = element_text(size = 12),
        panel.grid.major.y=element_blank()) +
  guides(color=guide_legend(reverse=TRUE))


EnvComp <- p + ggforce::facet_col(vars(Param), scales = 'free_y', space = 'free')


### Genetic component
m6 <- readRDS(file="outputs/models/P1/MOD6.rds")

# extract the posterior distribution of the parameters of interest
df <- m6 %>% broom::tidyMCMC(estimate.method = "median",conf.int = T,conf.level = 0.95) %>% 
  filter(str_detect(term, "^r_prov\\[|^r_mmQ1Q2Q3Q4Q5Q6\\[|^r_prov_clim")) %>% 
  mutate(Param = case_when(
    grepl("mmQ",term) ~ "Intercepts of the gene pools",
    grepl("r_prov\\[",term) ~ "Intercepts of the provenances",
    grepl("r_prov_clim",term) ~ "Intercepts of the climate-of-origin"))

# extract information on the main gene pool for each provenance
data <- readRDS(file="data/TrainP1.RDS")
data <- droplevels(data)
ps <- data %>%
  group_by(prov) %>% 
  summarise_at(vars(paste0(rep("Q",6),1:6)), mean) %>% 
  dplyr::rename(term=prov)
ps$GP <- colnames(ps[,2:7])[apply(ps[,2:7],1,which.max)]
ps$term <- as.factor(ps$term)

# Merge the two
prov <- df %>% 
  dplyr::filter(str_detect(term, "^r_prov")) %>% 
  mutate(term=case_when(startsWith(term, "r_prov[") ~str_sub(term,8,-12),
                        startsWith(term, "r_prov_clim[") ~ str_sub(term,13,-12))) %>% 
  left_join(ps[,c("term","GP")],by="term")

# >> Intercepts of the provenances
provint <- prov %>% filter(Param=="Intercepts of the provenances") %>% 
  group_by(term) %>%
  dplyr::mutate(meanperprov = mean(estimate))%>%
  ungroup()


pprov <- provint %>% 
  ggplot(aes(x = reorder(as.factor(term), meanperprov), y =estimate ,ymin = conf.low, ymax = conf.high,color=GP)) +
  coord_flip() +
  facet_grid(.~Param,scales="free_y", space = "free") +
  geom_pointinterval(position = position_dodge(width = .6),point_size=2.5,size=2.5) +
  ylab("") + xlab("") +
  theme_bw() +
  scale_color_manual(values=c("orangered3",
                             "gold2","darkorchid3",
                             "navyblue",
                             "turquoise2",
                             "green3")) +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=16),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=11),
        legend.position = "none",
        legend.background = element_rect(colour = "grey"),
        panel.grid.minor.y=element_blank(),
        strip.text.x = element_text(size = 12),
        panel.grid.major.y=element_blank())

# Intercepts of the climate-of-origin
provclim <- prov %>% filter(Param=="Intercepts of the climate-of-origin") %>% 
  group_by(term) %>%
  dplyr::mutate(meanperprov = mean(estimate))%>%
  ungroup()


pclim <- provclim %>% 
  ggplot(aes(x = reorder(as.factor(term), meanperprov), y =estimate ,ymin = conf.low, ymax = conf.high,color=GP)) +
  coord_flip() +
  facet_grid(.~Param,scales="free_y", space = "free") +
  #facet_grid(rows = vars(Param),scales="free_y", space = "free") +
  #facet_wrap(.~Param,scales="free_y",ncol=1) +
  geom_pointinterval(position = position_dodge(width = .6),point_size=2.5,size=2.5) +
  ylab("") + xlab("") +
  theme_bw() +
  scale_color_manual(values=c("orangered3",
                              "gold2","darkorchid3",
                              "navyblue",
                              "turquoise2",
                              "green3")) +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=16),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=11),
        legend.position = "none",
        legend.background = element_rect(colour = "grey"),
        panel.grid.minor.y=element_blank(),
        strip.text.x = element_text(size = 12),
        panel.grid.major.y=element_blank())



# Intercepts of the gene pools
gp <- df %>% dplyr::filter(str_detect(term, "mmQ")) %>% 
  mutate(term=str_sub(term,18,-12),
         GP=term,
         term=case_when(term=="Q1" ~ "NA",
                        term=="Q2" ~ "C",
                        term=="Q3" ~ "CS",
                        term=="Q4" ~ "FA",
                        term=="Q5" ~ "IA",
                        term=="Q6" ~ "SES")) %>% 
  mutate(GP=case_when(GP=="Q1" ~ "Northern Africa (NA)",
                      GP=="Q2" ~ "Corsica (C)",
                      GP=="Q3" ~ "Central Spain (CS)",
                      GP=="Q4" ~ "French Atlantic (FA)",
                      GP=="Q5" ~ "Iberian Atlantic (IA)",
                      GP=="Q6" ~ "South-eastern Spain (SES)")) %>% 
  
  group_by(term) %>%
  dplyr::mutate(meanperprov = mean(estimate))%>%
  ungroup()


gp$GP <- as.factor(gp$GP)
gp$GP <- factor(gp$GP,levels=c("French Atlantic (FA)",
                               "Iberian Atlantic (IA)",
                               "Corsica (C)",
                               "South-eastern Spain (SES)",
                               "Central Spain (CS)",
                               "Northern Africa (NA)"))

pgp <- gp %>% ggplot(aes(x = reorder(as.factor(term), meanperprov), 
                       y = estimate,ymin = conf.low, ymax = conf.high,color=GP)) +
  coord_flip() +
  facet_grid(.~Param,scales="free_y", space = "free") +
  #facet_grid(rows = vars(Param),scales="free_y", space = "free") +
  #facet_wrap(.~Param,scales="free_y",ncol=1) +
  geom_pointinterval(position = position_dodge(width = .6),point_size=2.5,size=2.5) +
  ylab("") + xlab("") +
  labs(color = "Gene pools") +
  theme_bw() +
  scale_color_manual(values=c("navyblue",
                              "turquoise2",
                              "gold2",
                              "green3",
                              "darkorchid3",
                              "orangered3")) +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=16),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=16),
        legend.position = "bottom",
        #legend.background = element_rect(colour = "grey"),
        panel.grid.minor.y=element_blank(),
        strip.text.x = element_text(size = 12),
        panel.grid.major.y=element_blank()) +
  guides(color=guide_legend(ncol=1,title.position="top", title.hjust = 0.5))

# Extract legend
legend <- get_legend(pgp)
pgp <- pgp + theme(legend.position = "none")

pgp <- plot_grid(pgp,legend,nrow=2)

GenComp <- plot_grid(pprov,pclim,pgp,nrow=1)


Fig3 <- plot_grid(plot_grid(VarPart,EnvComp,nrow=1),GenComp,nrow=2)

ggsave(Fig3, file="figs/manuscript/Fig3AfterReview.png",width = 16,height=16) 
ggsave(Fig3,file="figs/manuscript/Fig3AfterReview.pdf", width=16, height=16, labels = c("A) Variance partioning"),
       label_size = 20)

