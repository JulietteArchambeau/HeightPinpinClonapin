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
