########################################################################################################################"
#                                                                                                                      #
#                   Estimating gene-pool specific heritabilities with M5 model                                         #
#                                 Tables S27 and S28, Figure S13                                                       #
#                                                                                                                      #
#                                   Juliette Archambeau                                                                #
#                                       18/03/2022                                                                     #
#                                                                                                                      #
########################################################################################################################"

library(brms)       # CRAN v2.11.1
library(xtable)     # CRAN v1.8-4
library(tidyr)      # CRAN v1.0.0
library(readr)      # CRAN v1.3.1
library(ggpubr)     # CRAN v0.2.1
library(broom)      # CRAN v0.5.2
library(latex2exp)  # CRAN v0.4.0
library(stringr)    # CRAN v1.4.0
library(ggridges)   # CRAN v0.5.1
library(rethinking) # [github::rmcelreath/rethinking] v1.59
library(dplyr)      # CRAN v1.0.0
library(tibble)     # CRAN v2.1.3

mod <- readRDS(file="outputs/models/P1/MOD5.rds")

data <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv") %>%  dplyr::filter(P1=="train")

# Table S27   ####
# ========    "

# This table shows the outputs of a one-sided hypothesis testing on the probability 
# that the gene pool-specific total genetic variances in M5 are different.


clones <- paste0("clon",1:6)


for(i in clones[1:length(clones)]){
  for(j in clones[1:length(clones)]){
    if(i==j){next}
    hyp <- hypothesis(mod,paste0(i,"__Intercept<",j,"__Intercept"),class="sd")
    if(i=="clon1"&j=="clon2"){
      df <- as.data.frame(hyp$hypothesis)
    } else {
      df <- bind_rows(df,hyp$hypothesis)
    }
  }}



print(xtable(df, type = "latex",digits=3), 
      file = paste0("tables/Posteriors/M5_CompareTotGenVar.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})


# Figure S13 and Table S28   ####
# ========================   "

# a) Extracting the modeled residual variance:

sigma2 <- as.matrix(mod, pars = c("sigma"))^2
data.frame(Parameter="$\\sigma^{2}$",
           Median=median(sigma2),
           SD=sd(sigma2),
           LowCI=HPDI(sigma2,prob=0.95)[[1]],
           SudCI=HPDI(sigma2,prob=0.95)[[2]]) %>% 
  knitr::kable(digits = 4 ,format = "markdown")
var_res <- sigma2


# b) Extracting the gene pool-specific heritabilities

# >>> Northern-Africa gene pool
var_gen <- as.matrix(mod, pars = c("sd_clon1__Intercept"))^2
heritability <- var_gen / (var_gen + var_res)
df <- data.frame(Median=mean(heritability),
                 SD=sd(heritability),
                 InfCI=HPDI(heritability,prob=0.95)[[1]],
                 SupCI=HPDI(heritability,prob=0.95)[[2]])
row.names(df) <- "hNA"
hNA <- heritability[,1]

# >>> Corsican gene-pool
var_gen <- as.matrix(mod, pars = c("sd_clon2__Intercept"))^2
heritability <- var_gen / (var_gen + var_res)
df["hC",] <- c(mean(heritability),
               sd(heritability),
               HPDI(heritability,prob=0.95)[[1]],
               HPDI(heritability,prob=0.95)[[2]])
hC <- heritability[,1]

# >>> Central Spain gene pool
var_gen <- as.matrix(mod, pars = c("sd_clon3__Intercept"))^2
heritability <- var_gen / (var_gen + var_res)
df["hCS",] <- c(mean(heritability),
                sd(heritability),
                HPDI(heritability,prob=0.95)[[1]],
                HPDI(heritability,prob=0.95)[[2]])
hCS <- heritability[,1]

# >>> French Atlantic gene pool
var_gen <- as.matrix(mod, pars = c("sd_clon4__Intercept"))^2
heritability <- var_gen / (var_gen + var_res)
df["hFA",] <- c(mean(heritability),
                sd(heritability),
                HPDI(heritability,prob=0.95)[[1]],
                HPDI(heritability,prob=0.95)[[2]])
hFA <- heritability[,1]

# >>> Iberian Atlantic gene pool
var_gen <- as.matrix(mod, pars = c("sd_clon5__Intercept"))^2
heritability <- var_gen / (var_gen + var_res)
df["hIA",] <- c(mean(heritability),
                sd(heritability),
                HPDI(heritability,prob=0.95)[[1]],
                HPDI(heritability,prob=0.95)[[2]])
hIA <- heritability[,1]

# >>> South-eastern Spain gene pool
var_gen <- as.matrix(mod, pars = c("sd_clon6__Intercept"))^2
heritability <- var_gen / (var_gen + var_res)
df["hSES",] <- c(mean(heritability),
                 sd(heritability),
                 HPDI(heritability,prob=0.95)[[1]],
                 HPDI(heritability,prob=0.95)[[2]])
hSES <- heritability[,1]

# c) Creating and saving Table S28:
df <- df %>% mutate(Parameter = recode_factor(rownames(df),
                                              'hNA'='$H^{2}_{NA}$',
                                              'hC'='$H^{2}_{C}$',
                                              'hCS'='$H^{2}_{CS}$',
                                              'hFA'='$H^{2}_{FA}$',
                                              'hIA'='$H^{2}_{IA}$',
                                              'hSES'='$H^{2}_{SES}$')) %>%  
  remove_rownames() %>%
  dplyr::select(Parameter,Median,SD,InfCI,SupCI)


print(xtable(df, type = "latex",digits=3), 
      file = paste0("tables/Posteriors/M5_HeritabilityEstimates.tex"), 
      include.rownames=FALSE,
      sanitize.text.function = function(x) {x})


# d) Creating Figure S13:


POST <- data.frame(genepool=c(rep("hNA",6000),rep("hC",6000),rep("hCS",6000),rep("hFA",6000),rep("hIA",6000),rep("hSES",6000)),
                   value=c(hNA,hC,hCS,hFA,hIA,hSES))

POST <- POST %>%
  group_by(genepool) %>%
  dplyr::mutate(meanpergenepool = mean(value))%>%
  ungroup()

# >>>> Panel a
p <- ggplot()+
  geom_vline(xintercept = 0, col = "grey70") +
  geom_density_ridges(data  = POST, 
                      aes(x      = value,
                          y      = reorder(as.factor(genepool), meanpergenepool),
                          fill   = as.factor(genepool), vline_color = ..quantile..),
                      scale = 1, 
                      alpha = .6,
                      rel_min_height=c(.006),
                      size=0.8,
                      quantile_lines = TRUE, quantiles = c(0.025,0.5,0.975)) +
  
  scale_discrete_manual("vline_color",
                        values = c("blue", "red", "blue", "black"), 
                        breaks = c(1, 2),
                        labels = c("5th & 95th percentiles", "Median"),
                        name = NULL) +
  
  scale_fill_manual(values=c("gold2","darkorchid3","navyblue","turquoise2","orangered3","green3"), 
                    labels = c("Corsica (C)",
                               "Central Spain (CS)",
                               "French Atlantic (FA)",
                               "Iberian Atlantic (IA)",
                               "Northern Africa (NA)", 
                               "South-eastern Spain (SES)"),
                    name="Gene pools:") +
  
  scale_y_discrete(labels=c("hNA"=parse(text = TeX("$H_{NA}^{2}$")),
                            "hC"=parse(text = TeX("$H_{C}^{2}$")),
                            "hCS"=parse(text = TeX("$H_{CS}^{2}$")),
                            "hFA"=parse(text = TeX("$H_{FA}^{2}$")),
                            "hIA"=parse(text = TeX("$H_{IA}^{2}$")),
                            "hSES"=parse(text = TeX("$H_{SES}^{2}$")))) +
  labs(y  = "", 
       x  = TeX("Gene-pool specific broad-sense heritability estimates")) +
  
  theme_bw() + 
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=18), 
        legend.text = element_text(size=18,margin = margin(t = 10)),
        legend.title = element_text(size=20),
        legend.spacing.y = unit(1, 'cm')) 


# >>> Panel b

# >>>>>>> Calculating the partial breeding values for each gene pool

clon1 <- tidy(mod,parameters="r_clon1", intervals = FALSE)
clon1$clon <- str_sub(clon1$term,9,-12) 
clon1$term <- NULL
colnames(clon1) <- c("est1","std1","clon")

clon2 <- tidy(mod,parameters="r_clon2", intervals = FALSE)
clon2$clon <- str_sub(clon2$term,9,-12) 
clon2$term <- NULL
colnames(clon2) <- c("est2","std2","clon")

clon3 <- tidy(mod,parameters="r_clon3", intervals = FALSE)
clon3$clon <- str_sub(clon3$term,9,-12) 
clon3$term <- NULL
colnames(clon3) <- c("est3","std3","clon")

clon4 <- tidy(mod,parameters="r_clon4", intervals = FALSE)
clon4$clon <- str_sub(clon4$term,9,-12) 
clon4$term <- NULL
colnames(clon4) <- c("est4","std4","clon")

clon5 <- tidy(mod,parameters="r_clon5", intervals = FALSE)
clon5$clon <- str_sub(clon5$term,9,-12) 
clon5$term <- NULL
colnames(clon5) <- c("est5","std5","clon")

clon6 <- tidy(mod,parameters="r_clon6", intervals = FALSE)
clon6$clon <- str_sub(clon6$term,9,-12) 
clon6$term <- NULL
colnames(clon6) <- c("est6","std6","clon")

ps <- data %>% 
  dplyr::select(clon,Q1,Q2,Q3,Q4,Q5,Q6,max.Qvalue,max.Q) %>%  
  distinct()

g <- Reduce(
  function(x, y, ...) merge(x, y, ...), 
  list(clon1,clon2,clon3,clon4,clon5,clon6,ps)
)

g <- g %>% 
  dplyr::select(starts_with("est"),clon,max.Q) %>%
  pivot_longer(-c(clon,max.Q),names_to="PBV",values_to="est")

# >>>>>>>> Creating and saving Figure S13

pPBV <- ggplot(g, aes(x=max.Q, y=est, fill=PBV)) + 
  geom_boxplot(alpha=0.6) + theme_bw() +
  scale_x_discrete(labels=c("Q1"="Northern Africa",
                            "Q2"="Corsica",
                            "Q3"="Central Spain",
                            "Q4"="French Atlantic",
                            "Q5"="Iberian Atlantic",
                            "Q6"="South-eastern Spain")) + 
  scale_fill_manual(values=c("orangered3","gold2","darkorchid3","navyblue","turquoise2","green3"), 
                    labels = c("Northern Africa", 
                               "Corsica",
                               "Central Spain",
                               "French Atlantic",
                               "Iberian Atlantic",
                               "South-eastern Spain"),
                    name=TeX("Gene pools\ncontributing to the\npartial genetic values:")) + 
  labs(y        = "Means of the posterior distribution\nof the genetic values", 
       x        = TeX("Dominant gene pool for each genotype")) +
  theme(axis.text = element_text(size=14.5),
        axis.title = element_text(size=20), 
        legend.text = element_text(size=18),
        legend.title = element_text(size=20)) 


pPBVf <- pPBV + theme(legend.position = "none")

pp <-  ggarrange(p,pPBVf,labels=c("A","B"),nrow=2,heights = c(1, 1.3))

ggsave(pp,file="figs/SuppInfo/heritabilityM5.png",height=12,width=12)

