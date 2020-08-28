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

colnames(gp) <- c("CODE","NA","C","CS","FA","IA","SES")

data <- left_join(gp,prov)

write.csv(data,file="data/CoordProvGP.csv")
saveRDS(data,file="data/CoordProvGP.rds")
