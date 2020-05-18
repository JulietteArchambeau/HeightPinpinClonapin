###############################################################################################"
##################                                                      #######################"
##################           Script from Juliette Archambeau            #######################"
##################                    July  2019                        #######################"
##################               MERGING ALL DATABASES                  #######################"
##################            (pheno, climate and soil database)        #######################"
##################                                                      #######################"
###############################################################################################"

library(tidyverse)
library(dplyr)

pheno <- readRDS(file="data/PhenoData.RDS")

### Merge with site climate
clim.site <- readRDS(file="data/SiteMonthlyClimateData.RDS")

#pheno[,c(rownames(clim.site))] <- NA
for (i in rownames(clim.site)){
  pheno[pheno$site=="asturias"&pheno$age==10,i] <- clim.site[i,"ast2011dec"]
  pheno[pheno$site=="asturias"&pheno$age==21,i] <- clim.site[i,"ast2012nov"]
  pheno[pheno$site=="asturias"&pheno$age==37,i] <- clim.site[i,"ast2014mar"]
  
  pheno[pheno$site=="bordeaux"&pheno$age==13,i] <- clim.site[i,"bdx2012nov"]
  pheno[pheno$site=="bordeaux"&pheno$age==25,i] <- clim.site[i,"bdx2013nov"]
  pheno[pheno$site=="bordeaux"&pheno$age==37,i] <- clim.site[i,"bdx2014nov"]
  
  pheno[pheno$site=="caceres"&pheno$age==8,i] <- clim.site[i,"cac2011dec"]
  
  pheno[pheno$site=="madrid"&pheno$age==13,i] <- clim.site[i,"mad2011dec"]
  
  pheno[pheno$site=="portugal"&pheno$age==11,i] <- clim.site[i,"por2012jan"]
  pheno[pheno$site=="portugal"&pheno$age==15,i] <- clim.site[i,"por2012may"]
  pheno[pheno$site=="portugal"&pheno$age==20,i] <- clim.site[i,"por2012oct"]
  pheno[pheno$site=="portugal"&pheno$age==27,i] <- clim.site[i,"por2013may"]
}
sapply(pheno, function(x) sum(is.na(x)))


### Merge with prov climate 
clim.prov <- readRDS(file="data/ProvAnnualClimateData.RDS")
colnames(clim.prov) <- paste(colnames(clim.prov),"prov",sep="_") 
clim.prov <- as.tibble(cbind(rownames(clim.prov),data.frame(clim.prov,row.names = NULL))) %>% 
  rename(prov = "rownames(clim.prov)")
total <- merge(pheno, clim.prov, by="prov")

### Merge with prov and site soil 
site_soil <- readRDS(file="data/SiteSoilData.RDS")
site_soil <- as.tibble(cbind(rownames(site_soil),data.frame(site_soil,row.names = NULL))) %>% 
  rename(site = "rownames(site_soil)")
total <- merge(total, site_soil,by=c("site","latitude_site","longitude_site"))
prov_soil <- readRDS(file="data/ProvSoilData.RDS")
prov_soil <- as.tibble(cbind(rownames(prov_soil),data.frame(prov_soil,row.names = NULL))) %>% 
  rename(prov = "rownames(prov_soil)")
total <- merge(total, prov_soil, by=c("prov","latitude_prov","longitude_prov"))

#### Merge with SPEI
# prov_spei <- readRDS(file="data/prov_spei_data.RDS")
# prov_spei <- as.tibble(cbind(rownames(prov_spei),data.frame(prov_spei,row.names = NULL))) %>% 
#   rename(prov = "rownames(prov_spei)")
# total <- merge(total, prov_spei,by=c("prov"))



## Block as characters and not numeric
total$block <- as.character(total$block)
total$block_original <- as.character(total$block_original)





# Add population structure data ####

data <- read.delim("data/ClonapinBlups523IndPiMASSJuly2019.txt")
data <- data[,c("X","prov","Q1","Q2","Q3","Q4","Q5")]
data <- data %>% rename(clon=X)

# Checking NAs
sapply(data, function(x) sum(is.na(x)))

# Adding Q6
data$Q6 <- 1 - (data$Q1 + data$Q2 + data$Q3 + data$Q4 + data$Q5)

# Adding the cluster the most impt for each clone
for (i in 1:length(data$clon)){
  data[i,"max.Qvalue"] <- max(data[i,3:8])
  data[i,"max.Q"] <- names(data[i,3:8])[which.max(apply(data[i,3:8],MARGIN=2,max))]
}
data$max.Q <- as.factor(data$max.Q)


# 18 clones have no structure data
setdiff(unique(data[,"clon"]),unique(total[,"clon"]))
setdiff(unique(total[,"clon"]),unique(data[,"clon"]))
length(setdiff(unique(total[,"clon"]),unique(data[,"clon"])))

total <- merge(total,data,by=c("clon","prov"),all.x=TRUE)


#### SAVE ALL
saveRDS(total, file="data/AllDataPhenoClimSoil.RDS")

sapply(total, function(x) sum(is.na(x)))
