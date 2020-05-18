#################################################################################################"
####################                  February 2019                    ##########################"
####################     Computing the mean of climatic variables      ##########################"
###################             for the period 1901 - 2009             ##########################"
###################                in each provenance                  ##########################"
#################################################################################################"


# function to select dates included in a given time period
select.yr <- function(tab,period) c(3:116)[c(1901:2014)>=min(period) & c(1901:2014)<=max(period)]

# Which climatic data?
bioclim=paste("bio",c(1,2,5,6,12,13,14),sep="")
T.seas=paste0("tmean.",c("djf","mam","jja","son"))
P.seas=paste0("prec.",c("djf","mam","jja","son"))
pet=paste0("pet.",c("mean","min","max"))
ppet=paste0("ppet.",c("mean","min","max"))
eumedclim.vars=c(bioclim,T.seas,P.seas,pet,ppet)
coord <- readRDS("data/SiteProvCoord.RDS")

# keep only coordinates of the provenances
coord <- coord[c(1:35),]

df <- matrix(NA,dim(coord)[[1]],length(eumedclim.vars),dim=list(c(rownames(coord)),c(eumedclim.vars)))

for(k in eumedclim.vars) {
  nm.file=paste0("data/climate/extraction_",k,"_40points_eumedclim.RData")
  load(nm.file)
  clim.k <- clim.k[c(1:35),]
  df[,k]=apply(clim.k[,select.yr(tab=clim.k,1901:2009)],1,mean)
}
df <- as.data.frame(df)

saveRDS(df,file="data/ProvAnnualClimateData.RDS")
