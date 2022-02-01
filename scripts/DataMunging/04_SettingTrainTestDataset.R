###############################################################################################"
##################                                                      #######################"
##################                      July 2019                       #######################"
##################                TRAIN and TEST dataset                #######################"
##################                                                      #######################"
###############################################################################################"

library(tidyr)
library(dplyr)
library(compare)
library(stringr)

data <- readRDS(file="data/AllDataPhenoClimSoil.RDS")
data <- as.data.frame(data)

# Random terms as factors
randoms <- c("block","site","prov","clon","tree")
for (i in randoms){
  data[,i] <- as.factor( data[,i])
}
str(data[,c("block","site","prov","clon","tree")] )



#### Remove NAs in STRUCTURE data
#################################
# Comment: All the trees in the pop ROD have NAs for the structure... We are going to remove the entire provenance...
compare(data[data$prov=="ROD",],data[is.na(data$Q1)&data$prov=="ROD",])


# 12 populations have genotypes with NAs
length(unique(data$prov[is.na(data$Q1)]))
unique(data$prov[is.na(data$Q1)])

# 18 genotypes have NAs
length(unique(data$clon[is.na(data$Q1)]))


# removing NAs for population structure
sapply(data, function(x) sum(is.na(x)))
data <- data[!(is.na(data$Q1)),]
sapply(data, function(x) sum(is.na(x)))


#### Removing missing data for height
#####################################

data <- data[!(is.na(data$height)),]
sapply(data, function(x) sum(is.na(x)))


data <- as_tibble(data)


# !!!!  Important !!!!
# I forgot to use set.seed() and therefore the sampling of the different partitions is not replicable.
# Below are the code lines that have been run to generate the partitions.


#  P1 partition: sampling random observations
#############################################

sample <- sample.int(nrow(data), 0.75*nrow(data),replace=F)

datatest <- data[-sample,] # 25% of the data
datatrain <- data[sample,] # 75% of the data

#saveRDS(datatest, file="data/TestP1.RDS")
#saveRDS(datatrain, file="data/TrainP1.RDS")




#  P2 partition: sampling random provenances
############################################

selected.provs <- sample(unique(data$prov),28,replace = F)

train <- data[data$prov %in% selected.provs, ] # 28 provenances
test <- anti_join(data,train) # 6 provenances

#saveRDS(test, file="data/TestP2.RDS")
#saveRDS(train, file="data/TrainP2.RDS")



#  P3 partition: sampling provenances (not totally randomly) - 23/04/2020
#########################################################################

selected.provs <- c("PIA","TAM","ORI",as.character(sample(unique(data$prov),3,replace = F)))

test <- data[data$prov %in% selected.provs, ]
train <- anti_join(data,test)

#saveRDS(test, file="data/TestP3.RDS")
#saveRDS(train, file="data/TrainP3.RDS")



# We save data for the DRYAD repository associated with the AmNat paper:
########################################################################

# 1/ For more clarity, we keep only the variables used in the paper
# 2/ As we can not replicate the sampling of the three partitions (I forgot to use set.seed before partitioning the dataset), 
#    we add three variables (P1, P2, P3) indicating for each partition (P1, P2, P3) which observations belong to the test or train dataset.
# 3/ we save the data in comma-separated variable csv file as advised by AmNat: http://comments.amnat.org/2021/12/guidelines-for-archiving-code-with-data.html

dataDryad <- data %>% 
  dplyr::select(obs,tree,site,clon,prov,      # variables related to the experimental design
                latitude_site,longitude_site, # coord of the common gardens
                latitude_prov,longitude_prov, # coord of the provenances
                age,height,survival,
                pre_summer_min_site, pre_mean_1y_site, tmn_min_1y_site, tmx_max_1y_site, pre_max_1y_site, tmx_mean_1y_site, # climatic variables of the test sites
                bio1_prov, bio5_prov, bio12_prov, bio14_prov # climatic variables of the provenances
                ) %>% 
  dplyr::mutate(P1 = case_when(obs %in% (readRDS(file = "data/TrainP1.RDS") %>% dplyr::select(obs) %>% pull()) ~ "train",
                               obs %in% (readRDS(file = "data/TestP1.RDS") %>% dplyr::select(obs) %>% pull()) ~ "test"),
                P2 = case_when(obs %in% (readRDS(file = "data/TrainP2.RDS") %>% dplyr::select(obs) %>% pull()) ~ "train",
                               obs %in% (readRDS(file = "data/TestP2.RDS") %>% dplyr::select(obs) %>% pull()) ~ "test"),
                P3 = case_when(obs %in% (readRDS(file = "data/TrainP3.RDS") %>% dplyr::select(obs) %>% pull()) ~ "train",
                               obs %in% (readRDS(file = "data/TestP3.RDS") %>% dplyr::select(obs) %>% pull()) ~ "test"))

write.csv(dataDryad,
          file= paste0("data_DRYAD/HeightClimateSoilData_",nrow(dataDryad),"obs_",ncol(dataDryad),"variables.csv"),
          row.names=T)