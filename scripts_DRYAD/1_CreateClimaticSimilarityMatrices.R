##################################################################################################################
#                                                                                                                #
#        Creating covariance matrices describing the climatic similarity among provenances and test sites        #
#                                                                                                                #
#                                   Juliette Archambeau                                                          #
#                                       03/02/2022                                                               #
#                                                                                                                #
##################################################################################################################


# Detailed version of the script: https://github.com/JulietteArchambeau/HeightPinpinClonapin/blob/master/reports/ExplanatoryVariables/ClimSimMatrices.Rmd

# Packages:
library(readr)
library(matrixcalc)
library(Matrix)
library(tidyverse)


# Load data:
data <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv")



# I) Climatic similarity among test sites during the year preceding the measurements
# ==================================================================================

# create a column to identify each unique combination of site + year of measurement (= different age of the trees).
data$site_age <- paste0(data$site,data$age)

# selected variables:
  # `pre_summer_min_site` (*min.presummer* in the manuscript): minimum of the monthly precipitation during summer -June to September- in the test sites (°C)
  # `pre_mean_1y_site` (*mean.pre* in the manuscript): mean of the monthly precipitation in the test sites (mm)
  # `tmn_min_1y_site` (*min.tmn* in the manuscript):  minimum of the monthly minimum temperatures in the test sites (°C)
  # `tmx_max_1y_site` (*max.tmx* in the manuscript): maximum of the monthly maximum temperatures in the test sites (°C)
  # `pre_max_1y_site` (*max.pre* in the manuscript): maximum of the monthly precipitation in the test sites (mm)
  # `tmx_mean_1y_site` (*mean.tmax* in the manuscript): mean of monthly maximum temperatures in the test sites (°C)

select.var <- c("pre_summer_min_site","pre_mean_1y_site","tmn_min_1y_site",
                "tmx_max_1y_site","pre_max_1y_site","tmx_mean_1y_site")

# format data: one row per combination site/tree age and one column per climatic variable
dat <- unique(data[,c("site_age",select.var)])  # keep one row for each combination site/age
row.names(dat) <- dat$site_age                  # only climatic variables in column, and site/age as rownames
dat$site_age <- NULL 
dat <- scale(dat)                               # mean-center the selected variables

# calculate the variance-covariance matrix
varmat <-as.matrix(var(t(dat)))

# As the matrix has to be symetric and positive definite, we check it with the functions:
is.symmetric.matrix(varmat)                  # the matrix is symetric
is.positive.definite(varmat)                 # the matrix is not positive definite
varmat.pd <- nearPD(varmat)$mat              # to make the matrix positive definite
is.positive.definite(as.matrix(varmat.pd))   # the corrected matrix is now positive definite

write.csv(as.matrix(varmat.pd),
          file= paste0("data_DRYAD/VarCovMatSites.csv"),
          row.names=T)




# II) Climatic similarity among the 34 provenances (P1 partition)
# ===============================================================

# selected variables:
  # `bio1_prov` (*mean.temp* in the manuscript): the average of the annual daily mean temperature (°C).
  # `bio5_prov` (*max.temp* in the manuscript): the average of the maximum temperature of the warmest month (°C).
  # `bio12_prov` (*min.pre* in the manuscript): the average of the precipitation of the driest month (mm).
  # `bio14_prov` (*mean.pre* in the manuscript): the average of the annual precipitation (mm).

select.var <- paste0("bio",c(1,5,12,14),"_prov")

# format data: one row per provenance and one column per climatic variable
dat <- unique(data[,c("prov",select.var)])      # keep one row per provenance
row.names(dat) <- dat$prov                      # only climatic variables in column, and provenances as rownames
dat$prov <- NULL 
dat <- scale(dat)                               # mean-center the selected climatic variables


# calculate the variance-covariance matrix
varmat <-as.matrix(var(t(dat)))

# As the matrix has to be symetric and positive definite, we check it with the functions:
is.symmetric.matrix(varmat)                  # the matrix is symetric
is.positive.definite(varmat)                 # the matrix is not positive definite
varmat.pd <- nearPD(varmat)$mat              # to make the matrix positive definite
is.positive.definite(as.matrix(varmat.pd))   # the corrected matrix is now positive definite

write.csv(as.matrix(varmat.pd),
          file= paste0("data_DRYAD/VarCovMatProvenancesP1.csv"),
          row.names=T)

