############################################################################################"
##################  Script from Juliette Archambeau - July 2019      #######################"
##################         CLEANING PHENOTYPIC DATABASE              #######################"
############################################################################################"


# Libraries
library(readr)
library(tidyverse)
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)

# Import phenotypic data
data <- read_csv("data/CLONAPIN_all_data.csv")

# select the columns of interest
data <- data %>% select(region =Region, site, metap, prov, clon_original, block_original,
                        clon=clone, block, clones522, contains("_sur"),contains("_ht"),-contains("YEMA"))


# Remove 2 duplicated trees:
data[duplicated(data[,c("clon","block")]),c("clon","block")]
data[data$clon=="COC3"&data$block==37,]
data[data$clon_original=="STamr_F3_P4",]
data <- data[!(data$clon_original=="STamr_F3_P4"),]


data[data$clon=="TAM9"&data$block==10,c("BDX_htnov13","BDX_htnov14","BDX_htnov15")]
data[data$clon=="TAM9"&data$block==10&data$BDX_htnov13==60,
     c("BDX_htnov13","BDX_htnov14","BDX_htnov15")] <- colMeans(data[data$clon=="TAM9"&
                                                                      data$block==10,c("BDX_htnov13","BDX_htnov14","BDX_htnov15")])
data <- data[!(data$clon=="TAM9"&data$block==10&data$BDX_htnov13==54),]


### Replace prov "ASP" by "ARN" in Bordeaux because they are the same! Both codes refers to "Arenas de San Pedro" a population from Central Spain.
## Add 07/01/2019 -  see Marina de Miguel Vega's mail (07/01/2019)
# In the provenances:
data$prov[data$prov=="ASP"&data$site=="bordeaux"] <- "ARN"

# In the clones:
str_sub(data$clon[data$clon %like% "ASP"], 1,3) <- "ARN"


### A column to define a different IDs per tree (ID = clon_block)
data <- data %>% mutate(tree = NA) %>% unite(tree,clon,block, sep="_",remove=FALSE) 


# remane columns
data <- data %>% rename(AST_survdec11=AST_survdic11, CAC_survdec11=CAC_survdic11, MAD_survdec11=MAD_survdic11,
                        AST_htdec11=AST_htdic11_cm, AST_htnov12=AST_htnov12_cm, AST_htmar14=AST_htmar14_cm,
                        CAC_htdec11=CAC_htdic11,MAD_htdec11=MAD_htdic11)



### REMOVING ZOMBIE TREES
## Here, we delete measurements where trees were noted as "dead" when they were already noted as dead in previous measurements. 

### BORDEAUX
df <- subset(data, site=="bordeaux"&BDX_surv12==0&BDX_surv13==0&BDX_surv14==0&BDX_surv15==0)
data <- anti_join(data,df)
df$BDX_surv13 <- NA
df$BDX_surv14 <- NA
df$BDX_surv15 <- NA
data <- rbind(data,df)

df <- subset(data, site=="bordeaux"&BDX_surv13==0&BDX_surv14==0&BDX_surv15==0)
data <- anti_join(data,df)
df$BDX_surv14 <- NA
df$BDX_surv15 <- NA
data <- rbind(data,df)

df <- subset(data, site=="bordeaux"&BDX_surv14==0&BDX_surv15==0)
data <- anti_join(data,df)
df$BDX_surv15 <- NA
data <- rbind(data,df)

### ASTURIAS
df <- subset(data, site=="asturias"&AST_survdec11==0&AST_survnov12==0&AST_survmar14==0)
data <- anti_join(data,df)
df$AST_survnov12<- NA
df$AST_survmar14 <- NA
data <- rbind(data,df)

df <- subset(data, site=="asturias"&AST_survnov12==0&AST_survmar14==0)
data <- anti_join(data,df)
df$AST_survmar14 <- NA
data <- rbind(data,df)

### PORTUGAL
df <- subset(data, site=="portugal"&POR_survjan12==0&POR_survmay12==0&POR_survoct12==0&POR_survmay13==0)
data <- anti_join(data,df)
df$POR_survmay12<- NA
df$POR_survoct12 <- NA
df$POR_survmay13 <- NA
data <- rbind(data,df)

df <- subset(data, site=="portugal"&POR_survmay12==0&POR_survoct12==0&POR_survmay13==0)
data <- anti_join(data,df)
df$POR_survoct12 <- NA
df$POR_survmay13 <- NA
data <- rbind(data,df)

df <- subset(data, site=="portugal"&POR_survoct12==0&POR_survmay13==0)
data <- anti_join(data,df)
df$POR_survmay13 <- NA
data <- rbind(data,df)



## REMOVING THE "RESURRECTED TREES"
######## Some trees are noted as dead at a given measurement, but as alive in some of the following measurements. 
######## For these trees, we replace 0 (which means that the tree is dead) with 1 (the tree is alive) because we 
######## hypothesize that the trees were considered dead and might not be seen in the experiment (many tall grasses hiding the trees), 
######## but were actually alive (and noted as alive in the following measures). See Marina's email of 07/01/2019.

# BORDEAUX -> no trees concerned
df <- subset(data, site=="bordeaux"&(BDX_surv12==0|BDX_surv13==0|BDX_surv14==0)&BDX_surv15==1)
df <- subset(data, site=="bordeaux"&(BDX_surv12==0|BDX_surv13==0)&BDX_surv14==1)
df <- subset(data, site=="bordeaux"&BDX_surv12==0&BDX_surv13==1)

# ASTURIAS
df <- subset(data, site=="asturias"&(AST_survdec11==0|AST_survnov12==0)&AST_survmar14==1)
data <- anti_join(data,df)
df$AST_survdec11<- 1
df$AST_survnov12 <- 1
data <- rbind(data,df)

df <- subset(data, site=="asturias"&AST_survdec11==0&AST_survnov12==1)
data <- anti_join(data,df)
df$AST_survdec11<- 1
data <- rbind(data,df)

# PORTUGAL
df <- subset(data, site=="portugal"&(POR_survjan12==0|POR_survmay12==0|POR_survoct12==0)&POR_survmay13==1)
data <- anti_join(data,df)
df$POR_survjan12 <- 1
df$POR_survmay12<- 1
df$POR_survoct12 <- 1
data <- rbind(data,df)

df <- subset(data, site=="portugal"&(POR_survjan12==0|POR_survmay12==0)&POR_survoct12==1) # no trees concerned

df <- subset(data, site=="portugal"&POR_survjan12==0&POR_survmay12==1)
data <- anti_join(data,df)
df$POR_survjan12 <- 1
data <- rbind(data,df)


#*********** Be careful ! ****************************************************************************

# UNITS DIFFERENCES BETWEEN SITES
# In the file "CLONAPIN_all_data.csv", height measurements are:
# in cm for bordeaux and asturias
# in mm for caceres, madrid and portugal
# We are going to use mm for all sites:
col.cm <- c("AST_htdec11","AST_htnov12","AST_htmar14","BDX_htnov13","BDX_htnov14","BDX_htnov15")
data[,col.cm] <- data[,col.cm] * 10

# OUTLIERS 
# -> see script reports/00_outlier_pheno_data.Rmd
data[data$tree=="VAL16_7","AST_htmar14"] <- data[data$tree=="VAL16_7","AST_htmar14"] / 10
data[data$tree=="VER12_3","AST_htmar14"] <- data[data$tree=="VER12_3","AST_htmar14"] / 10

#******************************************************************************************************


# Trees were planted:
# 02/2011 in Asturias
# 04/2011 in Caceres
# 11/2011 in Bordeaux
# 11/2010 in Madrid
# 02/2011 in Portugal


###############################  ASTURIAS  ##################################################
####### december 2011
AST_dec11 <-  data[!(is.na(data$AST_survdec11)&is.na(data$AST_htdec11)),]
# 216 trees for which there was a value for survival but not for height. 
AST_dec11$survival <- AST_dec11$AST_survdec11
AST_dec11$height <- AST_dec11$AST_htdec11
AST_dec11$age <- 10
AST_dec11 <- AST_dec11 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)

####### november 2012
AST_nov12 <-  data[!(is.na(data$AST_survnov12)&is.na(data$AST_htnov12)),]
AST_nov12$survival <- AST_nov12$AST_survnov12
AST_nov12$height <- AST_nov12$AST_htnov12
AST_nov12$age <- 21
AST_nov12 <- AST_nov12 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)

####### march 2014
AST_mar14 <-  data[!(is.na(data$AST_survmar14)&is.na(data$AST_htmar14)),]
AST_mar14$survival <- AST_mar14$AST_survmar14
AST_mar14$height <- AST_mar14$AST_htmar14
AST_mar14$age <- 37
AST_mar14 <- AST_mar14 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)



###############################  BORDEAUX  ##################################################
####### november 2012
BDX_nov12 <-  data[!(is.na(data$BDX_surv12)),]
BDX_nov12$survival <- BDX_nov12$BDX_surv12
BDX_nov12$height <- NA
BDX_nov12$age <- 13
BDX_nov12 <- BDX_nov12 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)

####### november 2013
BDX_nov13 <-  data[!(is.na(data$BDX_surv13)&is.na(data$BDX_htnov13)),]
BDX_nov13$survival <- BDX_nov13$BDX_surv13
BDX_nov13$height <- BDX_nov13$BDX_htnov13
BDX_nov13$age <- 25
BDX_nov13 <- BDX_nov13 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)

####### november 2014
BDX_nov14 <-  data[!(is.na(data$BDX_surv14)&is.na(data$BDX_htnov14)),]
BDX_nov14$survival <- BDX_nov14$BDX_surv14
BDX_nov14$height <- BDX_nov14$BDX_htnov14
BDX_nov14$age <- 37
BDX_nov14 <- BDX_nov14 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)

####### november 2015 - but as climate data only goes until 2014, we can't include it them 
# BDX_nov15 <-  data[!(is.na(data$BDX_surv15)&is.na(data$BDX_htnov15)),]
# BDX_nov15$survival <- BDX_nov15$BDX_surv15
# BDX_nov15$height <- BDX_nov15$BDX_htnov15
# BDX_nov15$age <- 49
# BDX_nov15 <- BDX_nov15 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
#                                   clones522, age, survival, height)


###############################  CACERES  ##################################################
####### december 2011
CAC_dec11 <-  data[!(is.na(data$CAC_survdec11)&is.na(data$CAC_htdec11)),]
CAC_dec11$survival <- CAC_dec11$CAC_survdec11
CAC_dec11$height <- CAC_dec11$CAC_htdec11
CAC_dec11$age <- 8
CAC_dec11 <- CAC_dec11 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)


###############################  MADRID  ##################################################
####### december 2011
MAD_dec11 <-  data[!(is.na(data$MAD_survdec11)&is.na(data$MAD_htdec11)),]
MAD_dec11$survival <- MAD_dec11$MAD_survdec11
MAD_dec11$height <- MAD_dec11$MAD_htdec11
MAD_dec11$age <- 13
MAD_dec11 <- MAD_dec11 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)


###############################  PORTUGAL  ##################################################
####### january 2012
POR_jan12 <-  data[!(is.na(data$POR_survjan12)&is.na(data$POR_htjan12)),]
POR_jan12$survival <- POR_jan12$POR_survjan12
POR_jan12$height <- POR_jan12$POR_htjan12
POR_jan12$age <- 11
POR_jan12 <- POR_jan12 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)

####### may 2012
POR_may12 <-  data[!(is.na(data$POR_survmay12)&is.na(data$POR_htmay12)),]
POR_may12$survival <- POR_may12$POR_survmay12
POR_may12$height <- POR_may12$POR_htmay12
POR_may12$age <- 15
POR_may12 <- POR_may12 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)

####### october 2012
POR_oct12 <-  data[!(is.na(data$POR_survoct12)&is.na(data$POR_htoct12)),]
POR_oct12$survival <- POR_oct12$POR_survoct12
POR_oct12$height <- POR_oct12$POR_htoct12
POR_oct12$age <- 20
POR_oct12 <- POR_oct12 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)

####### may 2013
POR_may13 <-  data[!(is.na(data$POR_survmay13)&is.na(data$POR_htmay13)),]
POR_may13$survival <- POR_may13$POR_survmay13
POR_may13$height <- POR_may13$POR_htmay13
POR_may13$age <- 27
POR_may13 <- POR_may13 %>% select(region, site, metap, prov, clon_original, block_original, clon, tree, block,
                                  clones522, age, survival, height)


##########################################"
##########################################"

### Merge all
total <- rbind(AST_dec11,AST_nov12,AST_mar14,BDX_nov12,BDX_nov13,BDX_nov14,
               #BDX_nov15,
               CAC_dec11,MAD_dec11, POR_jan12,POR_may12,POR_oct12,POR_may13)

### quick check
table(total$age,total$site)

### A column to define a different IDs per obs (several osb per tree)
total <- total %>% mutate(obs = tree, num=NA) %>% group_by(tree) %>%  mutate(num = 1:n()) %>% 
  unite(obs, obs, num, sep="_") %>% select(1:5,7,10,6,9,8,obs,11:ncol(total))
length(unique(total$obs))

######################################"
## MISSING VALUES
## 8 trees have only NAs -> have not been included during the merging
names.total <- unique(total$tree) 
names.data <- unique(data$tree)
diff <- setdiff(names.data,names.total)

missing.values <- as.data.frame(matrix(NA,8,ncol(data),dimnames = list(c(diff),c(colnames(data)))))
for (i in diff){
  missing.values[i,] <- data[data$tree==i,]
}
######################################"


### Remove prov "SID" that is only in Bordeaux 
total <- total[!(total$prov=="SID"),]


### ADDING LATITUDE AND LONGITUDE
# Import & clean provenance coordinates
prov.coord <- read_csv("data/coordinates_provenances.csv")
prov.coord <- prov.coord[,c("CODE","ALTITUDE","LATITUDE","LONGITUDE")] %>%
  rename(prov=CODE,altitude_prov=ALTITUDE,latitude_prov=LATITUDE, longitude_prov=LONGITUDE) %>%
  slice(1:35) # remove the last line with only NAs (which remained from the original table)

# Merge provenance coordinates & phenotypic data
total <- merge(total,prov.coord,by="prov",all.x=T)
total <- as.tibble(total)

# Import & clean site coordinates
site.coord <- read_csv("data/coordinates_sites.csv")

# Merge site coordinates & phenotypic data
total <- merge(total,site.coord,by="site")
total <- as.tibble(total)


### Save ###################"
saveRDS(total,file="data/PhenoData.RDS")
############################"



