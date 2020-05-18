#################################################################################################"
####################                    July 2019                      ##########################"
####################     Computing the mean/min/max of climatic        ##########################"
###################         variables for given periods of time        ##########################"
###################                   for each site                    ##########################"
#################################################################################################"

library(SPEI)

## Loading files with sites coordinates
coord <- readRDS(file="data/SiteProvCoord.RDS")
coord <- coord[36:40,]

# function to select dates included in a given time period
select.yr <- function(tab,period) c(3:116)[c(1901:2014)>=min(period) & c(1901:2014)<=max(period)]


###### Extracting climate data for the period 2010 (January) - 2014 #####

## tmn = Min temperature
tmn=c(paste0("tmn",c('01','02','03','04','05','06','07','08','09','10','11','12')))
rep <- rep(tmn,5)
df.tmn <- matrix(NA,length(rep),5,dim=list(c(rep),c(rownames(coord[1:5,]))))
for(k in tmn) {
  nm.file=paste0("data/climate/extraction_",k,"_40points_eumedclim.RData")
  load(nm.file)
  df.tmn[rownames(df.tmn)==k,] <- t(subset(clim.k[c(36:40), select.yr(tab=clim.k,2010:2014)]))
}
df.tmn -> tmn
colnames(tmn) <- paste(colnames(tmn),"tmn",sep="_") 
rownames(tmn) <- 1:60

# tmx = Max temperature
tmx=c(paste0("tmx",c('01','02','03','04','05','06','07','08','09','10','11','12')))
rep <- rep(tmx,5)
df.tmx <- matrix(NA,length(rep),5,dim=list(c(rep),c(rownames(coord[1:5,]))))
for(k in tmx) {
  nm.file=paste0("data/climate/extraction_",k,"_40points_eumedclim.RData")
  load(nm.file)
  df.tmx[rownames(df.tmx)==k,] <- t(subset(clim.k[c(36:40), select.yr(tab=clim.k,2010:2014)]))
}
df.tmx -> tmx
colnames(tmx) <- paste(colnames(tmx),"tmx",sep="_") 
rownames(tmx) <-1:60

# pre = Precipitation
pre=c(paste0("pre",c('01','02','03','04','05','06','07','08','09','10','11','12')))
rep <- rep(pre,5)
df.pre <- matrix(NA,length(rep),5,dim=list(c(rep),c(rownames(coord[1:5,]))))
for(k in pre) {
  nm.file=paste0("data/climate/extraction_",k,"_40points_eumedclim.RData")
  load(nm.file)
  df.pre[rownames(df.pre)==k,] <- t(subset(clim.k[c(36:40), select.yr(tab=clim.k,2010:2014)]))
}
df.pre -> pre
colnames(pre) <- paste(colnames(pre),"pre",sep="_") 
rownames(pre) <- 1:60 

# merging tmn, tmx and pre
df <- cbind(pre,tmn,tmx)
rm(list=setdiff(ls(), c("df","coord")))



#### > adding 05/20/2019
#### Computing PET, the evapotranspiration potential using the latitude, the min and max temperatures and precipitation.
df <- as.data.frame(df)
df.pet <- sapply(row.names(coord), function(site){
  as.vector(SPEI::hargreaves(Tmin=df[,paste0(site,"_tmn")],Tmax=df[,paste0(site,"_tmx")],
                             lat=round(coord[site,"latitude"],4),Pre=df[,paste0(site,"_pre")]))
})
colnames(df.pet) <-  paste(colnames(df.pet),"pet", sep = "_")
df <- cbind(df,df.pet)




#### > adding 05/25/2019
#### Computing water balance: P-PET
df.ppet <- sapply(row.names(coord), function(site){
  df[,paste0(site,"_pre")] - df[,paste0(site,"_pet")]
})
colnames(df.ppet) <-  paste(colnames(df.ppet),"ppet", sep = "_")
df <- cbind(df,df.ppet)

saveRDS(df,file="data/TimeSeries.rds")

# print df
# months <- rep(c("Jan","Fev","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),5)
# vec <- c()
# for(i in 2010:2014){vec <- c(vec,rep(i,12))}
# row.names(df) <- paste(months,vec,sep="")
# round(df,2)
# df[,c(21:25)]

#### Computing mean, max and min of tmn, tmx and pre during the period between the planting date and
#### the measurement date. 

var.clim <- c("tmn","tmx","pre","pet","ppet") 
periods <- c("ast2011dec","ast2012nov","ast2014mar","bdx2012nov","bdx2013nov","bdx2014nov",
             "cac2011dec","mad2011dec","por2012jan","por2012may","por2012oct","por2013may")


# MEAN
## Mean monthly values over the period between the measurement date and one year before. 
clim.mean.1y <- matrix(NA,length(var.clim),length(periods),dim=list(var.clim,periods))

for (i in var.clim){
  clim.mean.1y[i,"ast2011dec"] <- mean(df[12:24,paste0("asturias_",i)])
  clim.mean.1y[i,"ast2012nov"] <- mean(df[23:35,paste0("asturias_",i)])
  clim.mean.1y[i,"ast2014mar"] <- mean(df[39:51,paste0("asturias_",i)])
  
  clim.mean.1y[i,"bdx2012nov"] <- mean(df[23:35,paste0("bordeaux_",i)])
  clim.mean.1y[i,"bdx2013nov"] <- mean(df[35:47,paste0("bordeaux_",i)])
  clim.mean.1y[i,"bdx2014nov"] <- mean(df[47:59,paste0("bordeaux_",i)])
  
  clim.mean.1y[i,"cac2011dec"] <- mean(df[12:24,paste0("caceres_",i)])
  
  clim.mean.1y[i,"mad2011dec"] <- mean(df[12:24,paste0("madrid_",i)])
  
  clim.mean.1y[i,"por2012jan"] <- mean(df[13:25,paste0("portugal_",i)])
  clim.mean.1y[i,"por2012may"] <- mean(df[17:29,paste0("portugal_",i)])
  clim.mean.1y[i,"por2012oct"] <- mean(df[12:34,paste0("portugal_",i)])
  clim.mean.1y[i,"por2013may"] <- mean(df[39:41,paste0("portugal_",i)])
}

rownames(clim.mean.1y) <- paste(rownames(clim.mean.1y),"mean_1y",sep="_") 


## Mean monthly values over the period between the measurement date and the previous measurement. - added 07/02/2019
clim.mean.int <- matrix(NA,length(var.clim),length(periods),dim=list(var.clim,periods))

for (i in var.clim){
  clim.mean.int[i,"ast2011dec"] <- mean(df[14:24,paste0("asturias_",i)])
  clim.mean.int[i,"ast2012nov"] <- mean(df[25:35,paste0("asturias_",i)])
  clim.mean.int[i,"ast2014mar"] <- mean(df[36:51,paste0("asturias_",i)])
  
  clim.mean.int[i,"bdx2012nov"] <- mean(df[22:35,paste0("bordeaux_",i)])
  clim.mean.int[i,"bdx2013nov"] <- mean(df[36:47,paste0("bordeaux_",i)])
  clim.mean.int[i,"bdx2014nov"] <- mean(df[48:59,paste0("bordeaux_",i)])
  
  clim.mean.int[i,"cac2011dec"] <- mean(df[16:24,paste0("caceres_",i)])
  
  clim.mean.int[i,"mad2011dec"] <- mean(df[11:24,paste0("madrid_",i)])
  
  clim.mean.int[i,"por2012jan"] <- mean(df[14:25,paste0("portugal_",i)])
  clim.mean.int[i,"por2012may"] <- mean(df[26:29,paste0("portugal_",i)])
  clim.mean.int[i,"por2012oct"] <- mean(df[30:34,paste0("portugal_",i)])
  clim.mean.int[i,"por2013may"] <- mean(df[35:41,paste0("portugal_",i)])
}

rownames(clim.mean.int) <- paste(rownames(clim.mean.int),"mean_int",sep="_") 




# MIN
## Minimum monthly values over the period between the measurement date and one year before. 
clim.min.1y <- matrix(NA,length(var.clim),length(periods),dim=list(var.clim,periods))

for (i in var.clim){
  clim.min.1y[i,"ast2011dec"] <- min(df[12:24,paste0("asturias_",i)])
  clim.min.1y[i,"ast2012nov"] <- min(df[23:35,paste0("asturias_",i)])
  clim.min.1y[i,"ast2014mar"] <- min(df[39:51,paste0("asturias_",i)])
  
  clim.min.1y[i,"bdx2012nov"] <- min(df[23:35,paste0("bordeaux_",i)])
  clim.min.1y[i,"bdx2013nov"] <- min(df[35:47,paste0("bordeaux_",i)])
  clim.min.1y[i,"bdx2014nov"] <- min(df[47:59,paste0("bordeaux_",i)])
  
  clim.min.1y[i,"cac2011dec"] <- min(df[12:24,paste0("caceres_",i)])
  
  clim.min.1y[i,"mad2011dec"] <- min(df[12:24,paste0("madrid_",i)])
  
  clim.min.1y[i,"por2012jan"] <- min(df[13:25,paste0("portugal_",i)])
  clim.min.1y[i,"por2012may"] <- min(df[17:29,paste0("portugal_",i)])
  clim.min.1y[i,"por2012oct"] <- min(df[12:34,paste0("portugal_",i)])
  clim.min.1y[i,"por2013may"] <- min(df[19:41,paste0("portugal_",i)])
}

rownames(clim.min.1y) <- paste(rownames(clim.min.1y),"min_1y",sep="_") 

## Minimum monthly values over the period between the measurement date and the previous measurement. - added 07/02/2019
clim.min.int <- matrix(NA,length(var.clim),length(periods),dim=list(var.clim,periods))

for (i in var.clim){
  clim.min.int[i,"ast2011dec"] <- min(df[14:24,paste0("asturias_",i)])
  clim.min.int[i,"ast2012nov"] <- min(df[25:35,paste0("asturias_",i)])
  clim.min.int[i,"ast2014mar"] <- min(df[36:51,paste0("asturias_",i)])
  
  clim.min.int[i,"bdx2012nov"] <- min(df[22:35,paste0("bordeaux_",i)])
  clim.min.int[i,"bdx2013nov"] <- min(df[36:47,paste0("bordeaux_",i)])
  clim.min.int[i,"bdx2014nov"] <- min(df[48:59,paste0("bordeaux_",i)])
  
  clim.min.int[i,"cac2011dec"] <- min(df[16:24,paste0("caceres_",i)])
  
  clim.min.int[i,"mad2011dec"] <- min(df[11:24,paste0("madrid_",i)])
  
  clim.min.int[i,"por2012jan"] <- min(df[14:25,paste0("portugal_",i)])
  clim.min.int[i,"por2012may"] <- min(df[26:29,paste0("portugal_",i)])
  clim.min.int[i,"por2012oct"] <- min(df[30:34,paste0("portugal_",i)])
  clim.min.int[i,"por2013may"] <- min(df[35:41,paste0("portugal_",i)])
}

rownames(clim.min.int) <- paste(rownames(clim.min.int),"min_int",sep="_") 



# MAX
## Maximum monthly values over the period between the measurement date and one year before. 
clim.max.1y <- matrix(NA,length(var.clim),length(periods),dim=list(var.clim,periods))

for (i in var.clim){
  clim.max.1y[i,"ast2011dec"] <- max(df[12:24,paste0("asturias_",i)])
  clim.max.1y[i,"ast2012nov"] <- max(df[23:35,paste0("asturias_",i)])
  clim.max.1y[i,"ast2014mar"] <- max(df[39:51,paste0("asturias_",i)])
  
  clim.max.1y[i,"bdx2012nov"] <- max(df[23:35,paste0("bordeaux_",i)])
  clim.max.1y[i,"bdx2013nov"] <- max(df[35:47,paste0("bordeaux_",i)])
  clim.max.1y[i,"bdx2014nov"] <- max(df[47:59,paste0("bordeaux_",i)])
  
  clim.max.1y[i,"cac2011dec"] <- max(df[12:24,paste0("caceres_",i)])
  
  clim.max.1y[i,"mad2011dec"] <- max(df[12:24,paste0("madrid_",i)])
  
  clim.max.1y[i,"por2012jan"] <- max(df[13:25,paste0("portugal_",i)])
  clim.max.1y[i,"por2012may"] <- max(df[17:29,paste0("portugal_",i)])
  clim.max.1y[i,"por2012oct"] <- max(df[12:34,paste0("portugal_",i)])
  clim.max.1y[i,"por2013may"] <- max(df[19:41,paste0("portugal_",i)])
}

rownames(clim.max.1y) <- paste(rownames(clim.max.1y),"max_1y",sep="_") 

## Maximum monthly values over the period between the measurement date and the previous measurement. - added 07/02/2019
clim.max.int <- matrix(NA,length(var.clim),length(periods),dim=list(var.clim,periods))

for (i in var.clim){
  clim.max.int[i,"ast2011dec"] <- max(df[14:24,paste0("asturias_",i)])
  clim.max.int[i,"ast2012nov"] <- max(df[25:35,paste0("asturias_",i)])
  clim.max.int[i,"ast2014mar"] <- max(df[36:51,paste0("asturias_",i)])
  
  clim.max.int[i,"bdx2012nov"] <- max(df[22:35,paste0("bordeaux_",i)])
  clim.max.int[i,"bdx2013nov"] <- max(df[36:47,paste0("bordeaux_",i)])
  clim.max.int[i,"bdx2014nov"] <- max(df[48:59,paste0("bordeaux_",i)])
  
  clim.max.int[i,"cac2011dec"] <- max(df[16:24,paste0("caceres_",i)])
  
  clim.max.int[i,"mad2011dec"] <- max(df[11:24,paste0("madrid_",i)])
  
  clim.max.int[i,"por2012jan"] <- max(df[14:25,paste0("portugal_",i)])
  clim.max.int[i,"por2012may"] <- max(df[26:29,paste0("portugal_",i)])
  clim.max.int[i,"por2012oct"] <- max(df[30:34,paste0("portugal_",i)])
  clim.max.int[i,"por2013may"] <- max(df[35:41,paste0("portugal_",i)])
}

rownames(clim.max.int) <- paste(rownames(clim.max.int),"max_int",sep="_") 



# Add a column with minimum precipitation in summer (june -> september)  -- modif done 05/29/2019
# variables describing the summer drought of the previous year
clim.pre.summer <- matrix(NA,1,length(periods),dim=list("pre",periods))

for (i in "pre"){
  clim.pre.summer[i,"ast2011dec"] <- min(df[18:21,paste0("asturias_",i)])
  clim.pre.summer[i,"ast2012nov"] <- min(df[30:33,paste0("asturias_",i)])
  clim.pre.summer[i,"ast2014mar"] <- min(df[42:45,paste0("asturias_",i)])
  
  clim.pre.summer[i,"bdx2012nov"] <- min(df[30:33,paste0("bordeaux_",i)])
  clim.pre.summer[i,"bdx2013nov"] <- min(df[42:45,paste0("bordeaux_",i)])
  clim.pre.summer[i,"bdx2014nov"] <- min(df[54:57,paste0("bordeaux_",i)])
  
  clim.pre.summer[i,"cac2011dec"] <- min(df[18:21,paste0("caceres_",i)])
  
  clim.pre.summer[i,"mad2011dec"] <- min(df[18:21,paste0("madrid_",i)])
  
  clim.pre.summer[i,"por2012jan"] <- min(df[18:21,paste0("portugal_",i)])
  clim.pre.summer[i,"por2012may"] <- min(df[18:21,paste0("portugal_",i)])
  clim.pre.summer[i,"por2012oct"] <- min(df[30:33,paste0("portugal_",i)])
  clim.pre.summer[i,"por2013may"] <- min(df[30:33,paste0("portugal_",i)])
}
rownames(clim.pre.summer) <- paste(rownames(clim.pre.summer),"_summer_min",sep="") 


total <- rbind(clim.mean.1y,clim.min.1y,clim.max.1y,clim.mean.int,clim.min.int,clim.max.int,clim.pre.summer)
total <- as.data.frame(total)
rownames(total) <- paste(rownames(total),"site",sep="_") 

#################################################"
saveRDS(total, file="data/SiteMonthlyClimateData.RDS")
#################################################"
