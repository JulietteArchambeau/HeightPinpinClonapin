########################################################################################################################"
#                                                                                                                      #
#                 Preparing the test datatests of the P1, P2 and P3 partitions                                         #
#                                                                                                                      #
#                                   Juliette Archambeau                                                                #
#                                       10/03/2022                                                                     #
#                                                                                                                      #
########################################################################################################################"


library(dplyr)

# Important: In the test datasets, the explanatory variables have to be normalized with the mean and variance of the *train* datasets
# https://stackoverflow.com/questions/49444262/normalize-data-before-or-after-split-of-training-and-testing-data
# https://sebastianraschka.com/faq/docs/scale-training-test.html


# P1 (sampling random observations. Test data set of 25% observations)  ####
############################################################################


test <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv") %>%  dplyr::filter(P1=="test")
train <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv") %>%  dplyr::filter(P1=="train")

# AGE
test$age.sc <- (test$age - mean(train$age)) / sd(train$age)


# Population structure
colnames(test)[colnames(test) %in% c(paste0("Q",rep(1:6)))] <- c(paste0("prop_Q",rep(1:6)))
head(test[,c(paste0("prop_Q",rep(1:6)))])
test$Q1 <- "Q1"
test$Q2 <- "Q2"
test$Q3 <- "Q3"
test$Q4 <- "Q4"
test$Q5 <- "Q5"
test$Q6 <- "Q6"
sum(test[,c(paste0("prop_Q",rep(1:6)))]<0)
filter_at(test,c(paste0("prop_Q",rep(1:6))),any_vars(. < 0))
test$prop_Q6[test$prop_Q6<0] <- 0


# Intercepts of the 6 weighed Genomic relationship matrices
test$clon1 <- test$clon
test$clon2 <- test$clon
test$clon3 <- test$clon
test$clon4 <- test$clon
test$clon5 <- test$clon
test$clon6 <- test$clon


## Intercepts of the site climatic similarity
test$site_age <- paste0(test$site,test$age)


# Intercepts of the provenance climatic similarity
test$prov_clim <- test$prov


# Count of PEAs
snp.counts <- readRDS(file="data/CountPEAs.RDS")
test <- merge(test,snp.counts,by="clon")

# Region-specific PEAs
test$rPEA <- NA 
test$rPEA[test$site=="asturias"|test$site=="portugal"] <- test$count_ibatl_350[test$site=="asturias"|test$site=="portugal"]
test$rPEA[test$site=="caceres"|test$site=="madrid"] <- test$count_med_350[test$site=="madrid"|test$site=="caceres"]
test$rPEA[test$site=="bordeaux"] <- test$count_fratl_350[test$site=="bordeaux"]

train <- merge(train,snp.counts,by="clon")
train$rPEA <- NA 
train$rPEA[train$site=="asturias"|train$site=="portugal"] <- train$count_ibatl_350[train$site=="asturias"|train$site=="portugal"]
train$rPEA[train$site=="caceres"|train$site=="madrid"] <- train$count_med_350[train$site=="madrid"|train$site=="caceres"]
train$rPEA[train$site=="bordeaux"] <- train$count_fratl_350[train$site=="bordeaux"]

rm(snp.counts)

test$rPEA.sc <-  (test$rPEA - mean(train$rPEA)) / sd(train$rPEA)

# Global PEAs
test$gPEA.sc <- (test$count_all_350 - mean(train$count_all_350)) / sd(train$count_all_350)


# Site climatic variables
test$pre_mean_1y_site.sc <- (test$pre_mean_1y_site - mean(train$pre_mean_1y_site)) / sd(train$pre_mean_1y_site)
test$pre_summer_min_site.sc <- (test$pre_summer_min_site - mean(train$pre_summer_min_site)) / sd(train$pre_summer_min_site)
test$tmn_min_1y_site.sc <- (test$tmn_min_1y_site - mean(train$tmn_min_1y_site)) / sd(train$tmn_min_1y_site)
test$tmx_max_1y_site.sc <- (test$tmx_max_1y_site - mean(train$tmx_max_1y_site)) / sd(train$tmx_max_1y_site)
test$tmx_mean_1y_site.sc <- (test$tmx_mean_1y_site - mean(train$tmx_mean_1y_site)) / sd(train$tmx_mean_1y_site)


# Provenance climatic variables
test$bio14_prov.sc <- (test$bio14_prov - mean(train$bio14_prov)) / sd(train$bio14_prov)
test$bio5_prov.sc <- (test$bio5_prov - mean(train$bio5_prov)) / sd(train$bio5_prov)

write.csv(test, file= paste0("data_DRYAD/TestP1prepared.csv"), row.names = F)



# P2 (randomly selected provenances - 6 provenances in the test data set)   ####
################################################################################

test <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv") %>%  dplyr::filter(P2=="test")
train <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv") %>%  dplyr::filter(P2=="train")

# AGE
test$age.sc <- (test$age - mean(train$age)) / sd(train$age)


# Population structure
colnames(test)[colnames(test) %in% c(paste0("Q",rep(1:6)))] <- c(paste0("prop_Q",rep(1:6)))
head(test[,c(paste0("prop_Q",rep(1:6)))])
test$Q1 <- "Q1"
test$Q2 <- "Q2"
test$Q3 <- "Q3"
test$Q4 <- "Q4"
test$Q5 <- "Q5"
test$Q6 <- "Q6"
sum(test[,c(paste0("prop_Q",rep(1:6)))]<0)
filter_at(test,c(paste0("prop_Q",rep(1:6))),any_vars(. < 0))
test$prop_Q6[test$prop_Q6<0] <- 0


# Intercepts of the 6 weighed Genomic relationship matrices
test$clon1 <- test$clon
test$clon2 <- test$clon
test$clon3 <- test$clon
test$clon4 <- test$clon
test$clon5 <- test$clon
test$clon6 <- test$clon


## Intercepts of the site climatic similarity
test$site_age <- paste0(test$site,test$age)


# Intercepts of the provenance climatic similarity
test$prov_clim <- test$prov


# Count of the SNPs
snp.counts <- readRDS(file="data/CountPEAs.RDS")
test <- merge(test,snp.counts,by="clon")

# Region-specific PEAs
test$rPEA <- NA 
test$rPEA[test$site=="asturias"|test$site=="portugal"] <- test$count_ibatl_350[test$site=="asturias"|test$site=="portugal"]
test$rPEA[test$site=="caceres"|test$site=="madrid"] <- test$count_med_350[test$site=="madrid"|test$site=="caceres"]
test$rPEA[test$site=="bordeaux"] <- test$count_fratl_350[test$site=="bordeaux"]

train <- merge(train,snp.counts,by="clon")
train$rPEA <- NA 
train$rPEA[train$site=="asturias"|train$site=="portugal"] <- train$count_ibatl_350[train$site=="asturias"|train$site=="portugal"]
train$rPEA[train$site=="caceres"|train$site=="madrid"] <- train$count_med_350[train$site=="madrid"|train$site=="caceres"]
train$rPEA[train$site=="bordeaux"] <- train$count_fratl_350[train$site=="bordeaux"]

rm(snp.counts)

test$rPEA.sc <-  (test$rPEA - mean(train$rPEA)) / sd(train$rPEA)

# Global PEAs
test$gPEA.sc <- (test$count_all_350 - mean(train$count_all_350)) / sd(train$count_all_350)


# Site climatic variables
test$pre_mean_1y_site.sc <- (test$pre_mean_1y_site - mean(train$pre_mean_1y_site)) / sd(train$pre_mean_1y_site)
test$pre_summer_min_site.sc <- (test$pre_summer_min_site - mean(train$pre_summer_min_site)) / sd(train$pre_summer_min_site)
test$tmn_min_1y_site.sc <- (test$tmn_min_1y_site - mean(train$tmn_min_1y_site)) / sd(train$tmn_min_1y_site)
test$tmx_max_1y_site.sc <- (test$tmx_max_1y_site - mean(train$tmx_max_1y_site)) / sd(train$tmx_max_1y_site)
test$tmx_mean_1y_site.sc <- (test$tmx_mean_1y_site - mean(train$tmx_mean_1y_site)) / sd(train$tmx_mean_1y_site)


# Provenance climatic variables
test$bio14_prov.sc <- (test$bio14_prov - mean(train$bio14_prov)) / sd(train$bio14_prov)
test$bio5_prov.sc <- (test$bio5_prov - mean(train$bio5_prov)) / sd(train$bio5_prov)

write.csv(test, file= paste0("data_DRYAD/TestP2prepared.csv"), row.names = F)





# P3 (selected provenances - 6 provenances in the test data set)     ####
#########################################################################

test <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv") %>%  dplyr::filter(P3=="test")
train <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv") %>%  dplyr::filter(P3=="train")

# AGE
test$age.sc <- (test$age - mean(train$age)) / sd(train$age)


# Population structure
colnames(test)[colnames(test) %in% c(paste0("Q",rep(1:6)))] <- c(paste0("prop_Q",rep(1:6)))
head(test[,c(paste0("prop_Q",rep(1:6)))])
test$Q1 <- "Q1"
test$Q2 <- "Q2"
test$Q3 <- "Q3"
test$Q4 <- "Q4"
test$Q5 <- "Q5"
test$Q6 <- "Q6"
sum(test[,c(paste0("prop_Q",rep(1:6)))]<0)
filter_at(test,c(paste0("prop_Q",rep(1:6))),any_vars(. < 0))
test$prop_Q6[test$prop_Q6<0] <- 0


# Intercepts of the 6 weighed Genomic relationship matrices
test$clon1 <- test$clon
test$clon2 <- test$clon
test$clon3 <- test$clon
test$clon4 <- test$clon
test$clon5 <- test$clon
test$clon6 <- test$clon


## Intercepts of the site climatic similarity
test$site_age <- paste0(test$site,test$age)


# Intercepts of the provenance climatic similarity
test$prov_clim <- test$prov


# Count of the SNPs
snp.counts <- readRDS(file="data/CountPEAs.RDS")
test <- merge(test,snp.counts,by="clon")


# Region-specific PEAs
test$rPEA <- NA 
test$rPEA[test$site=="asturias"|test$site=="portugal"] <- test$count_ibatl_350[test$site=="asturias"|test$site=="portugal"]
test$rPEA[test$site=="caceres"|test$site=="madrid"] <- test$count_med_350[test$site=="madrid"|test$site=="caceres"]
test$rPEA[test$site=="bordeaux"] <- test$count_fratl_350[test$site=="bordeaux"]

train <- merge(train,snp.counts,by="clon")
train$rPEA <- NA 
train$rPEA[train$site=="asturias"|train$site=="portugal"] <- train$count_ibatl_350[train$site=="asturias"|train$site=="portugal"]
train$rPEA[train$site=="caceres"|train$site=="madrid"] <- train$count_med_350[train$site=="madrid"|train$site=="caceres"]
train$rPEA[train$site=="bordeaux"] <- train$count_fratl_350[train$site=="bordeaux"]

rm(snp.counts)

test$rPEA.sc <-  (test$rPEA - mean(train$rPEA)) / sd(train$rPEA)


# Global PEAs
test$gPEA.sc <- (test$count_all_350 - mean(train$count_all_350)) / sd(train$count_all_350)


# Site climatic variables
test$pre_mean_1y_site.sc <- (test$pre_mean_1y_site - mean(train$pre_mean_1y_site)) / sd(train$pre_mean_1y_site)
test$pre_summer_min_site.sc <- (test$pre_summer_min_site - mean(train$pre_summer_min_site)) / sd(train$pre_summer_min_site)
test$tmn_min_1y_site.sc <- (test$tmn_min_1y_site - mean(train$tmn_min_1y_site)) / sd(train$tmn_min_1y_site)
test$tmx_max_1y_site.sc <- (test$tmx_max_1y_site - mean(train$tmx_max_1y_site)) / sd(train$tmx_max_1y_site)
test$tmx_mean_1y_site.sc <- (test$tmx_mean_1y_site - mean(train$tmx_mean_1y_site)) / sd(train$tmx_mean_1y_site)

# Provenance climatic variables
test$bio14_prov.sc <- (test$bio14_prov - mean(train$bio14_prov)) / sd(train$bio14_prov)
test$bio5_prov.sc <- (test$bio5_prov - mean(train$bio5_prov)) / sd(train$bio5_prov)


# saveRDS(test, file="data/TestP3prepared.RDS")
write.csv(test, file= paste0("data_DRYAD/TestP3prepared.csv"), row.names = F)




