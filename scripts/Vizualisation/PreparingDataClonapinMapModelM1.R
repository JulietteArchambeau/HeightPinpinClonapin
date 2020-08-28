# Table for Qgis

library(dplyr)

prov <- read.csv(file="data/coordinates_provenances.csv")

gp <- readRDS(file="data/AllDataPhenoClimSoil.RDS")

# removing missing structure data
gp <- gp[!(is.na(gp$Q1)),]

# removing missing height data
gp <- gp[!(is.na(gp$height)),]

gp <- droplevels(gp) %>%
  group_by(prov) %>% 
  summarise_at(vars(paste0(rep("Q",6),1:6)), mean)  # %>% mutate_at(2:7,funs(round(., 3))) 

for (i in 1:length(gp$prov)){
  gp[i,"mainGP"] <- names(gp[i,2:7])[which.max(apply(gp[i,2:7],MARGIN=2,max))]
}

colnames(gp) <- c("CODE","NA","C","CS","FA","IA","SES","mainGP")

data <- left_join(gp,prov)

write.csv(data,file="data/CoordProvGPbis.csv")
