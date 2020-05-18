##################################################################################################"
####################  Script from Thibault Fréjaville - November 2018   ##########################"
#################### Monthly climate data from geo coordinates (WGS 84) ##########################"
##################################################################################################"


require(rgdal)
require(raster)
require(ncdf4)
require(parallel)
require(dplyr)
library(readr)

#########################################"
### function to extract climate data  ###
#########################################"

##### ARGUMENTS
### 'xy', matrix or data.frame (mandatory). Sampling of geographical coordinates to extract climate data 
#   (two-columns coordinates of longitude and latitude in WGS84).
### 'clim.var', character (optional). Vector containing the selected EuMedClim climatic variables. 
#   If NULL (by default), all the available variables are selected.
### 'radius', numeric (optional). Radius (in meters) of a buffer around each geographical coordinate. 
#   Pixel values within the buffer are averaged. If the distance between the geographical coordinate and the 
#   center of the pixel is less than or equal to dimension of the buffer, the pixel is included. 
#   If NULL (by default), no buffer is applied.
### 'dir.in', character (mandatory). Full directory path to the folder where climate files are stored.
#   If files are missing, set download=T to download them from http://gentree.data.inra.fr/climate/.
### 'download', logical (optional). If TRUE, required climate files will be downloaded 
#   from http://gentree.data.inra.fr/climate/ and stored (full path of the destination folder given in 'dir.in'). 
#   FALSE by default.
### 'dir.out', character (optional). Full directory path to the folder to save climate data (RData format).
#   If more than one climatic variable is selected, one file will be created for each. If NULL (by default), 
#   the extracted climate data is not saved.
### 'n.core', numeric (optional). Number of CPU cores to use (no parallelization by default; not available 
#   on Windows). Useful when selecting several climatic variables. Use 'detectCores()-1' to know the number 
#   of available cores.

get.climate<-function(xy,clim.var="monthly",period=1901:2014,fun=NULL,radius=NULL,dir.in,download=T,dir.out,n.core) {
  
  # monthly variables provided by EuMedClim
  monthly.data=c() ; for(md in c('tmn','tmx','pre')) monthly.data=c(monthly.data,paste0(md,c('01','02','03','04','05','06','07','08','09','10','11','12')))
  
  # verify conditions
  if(ncol(xy)!=2) stop('you should provide a two-column matrix of longitude and latitude data')
  if(!is.data.frame(xy)) xy=as.data.frame(xy)
  
  if(is.null(clim.var) | clim.var=='monthly') clim.var=monthly.data else
    for(ki in clim.var) if(all(ki!=eumedclim.vars))
      stop(paste('the climatic variable',ki,'is not provided by EuMedClim'))
  
  if(is.null(period))
    stop("the argument 'period' is empty. Select years within 1901-2014") else
      if(min(period)<1901 | max(period)>2014)
        stop("error in the 'period' argument. Select years within 1901-2014")
  
  if(length(period)>1 & !is.null(fun) & all(fun!=c('mean','min','max')))
    stop("error in the 'fun' argument.")
  
  if(is.null(dir.in))
    stop("the argument 'dir.in' is empty. Precise the full path of climate data folder (or where to store downloaded files)") else
      if(substr(dir.in,nchar(dir.in),nchar(dir.in))!='/')
        dir.in=paste0(dir.in,'/')
  
  if(is.null(dir.out))
    warning("the argument 'dir.out' is empty. To save climate data, precise the full path of destination") else
      if(substr(dir.out,nchar(dir.out),nchar(dir.out))!='/')
        dir.out=paste0(dir.out,'/')
  
  # minimum radius
  if(!is.null(radius)) if(radius<5000) radius=5000
  
  # tiles
  tiles=list(c(-20,0,20,45),c(0,20,20,45),c(20,40,20,45),c(40,60,20,45),
             c(-20,0,45,72),c(0,20,45,72),c(20,40,45,72),c(40,60,45,72))
  tiles=lapply(tiles,extent)
  # tiles to be selected
  tiles.in=c()
  for(t in 1:length(tiles))
    tiles.in[t] = any(xy[,1]>=xmin(tiles[[t]]) & xy[,1]<xmax(tiles[[t]]) & xy[,2]>=ymin(tiles[[t]]) & xy[,2]<ymax(tiles[[t]]))
  
  # compute list of time series dataframes (one dataframe by parameter with the geographical coordinates in rows and yearly climate values in columns)
  list.ts=mclapply(as.list(clim.var),function(k) {
    
    clim.k=matrix(NA,nrow(xy),114,dimnames=list(rownames(xy),1901:2014))
    
    for(t in c(1:8)[tiles.in]) {
      
      # tile t
      tile=tiles[[t]]
      
      # corresponding geographical coordinates
      in.tile <- xy[,1]>=xmin(tile) & xy[,1]<xmax(tile) & xy[,2]>=ymin(tile) & xy[,2]<ymax(tile)
      
      # get file
      file.nm=paste0(k,"_1901-2014_1km_lon_",xmin(tile),"_",xmax(tile),"_lat_",ymin(tile),"_",ymax(tile),"_eumedclim.tif")
      
      # find directory folder for variable k
      dir.in.k=dir.in
      if(all(list.files(dir.in)!=file.nm))
        if(any(list.dirs(dir.in,full.names=F)==k)) dir.in.k=paste0(dir.in,k,'/')
      
      # download file
      if(all(list.files(dir.in.k)!=file.nm))
        
        if(download)
          
          download.file(url=paste0('http://gentree.data.inra.fr/climate/datasets/',k,'/',file.nm),destfile=paste0(dir.in.k,file.nm),method='auto',cacheOK=F) else
            
            stop(paste("the following file",file.nm,"was not found in 'dir.in'. Verify the 'dir.in' argument or set download=T to automatically download climate files"))
      
      # load file
      stk=stack(paste0(dir.in.k,file.nm))
      names(stk)=1901:2014
      
      # extract yearly data for each geographical coordinate by applying a buffer (if 'radius' argument is not NULL, it returns the average value)
      clim.k.t=raster::extract(stk,xy[in.tile,],buffer=radius,fun=mean)
      clim.k[in.tile,]=clim.k.t
      
      # release memory
      rm(stk,clim.k.t)
      
    }
    
    # convert units (to °C or mm)
    clim.k=0.1*round(clim.k)
    
    # add coordinates
    clim.k=data.frame(xy,clim.k)
    
    # save
    if(!is.null(dir.out) & is.null(fun)) {
      
      nm.data=paste0(dir.out,"extraction_",k,"_",nrow(xy),"points",ifelse(!is.null(radius),paste0("_buffer",radius),""),"_eumedclim.RData")
      save(clim.k,file=nm.data)
      
    }
    
    return(clim.k)
    
  },mc.cores=n.core)
  
  # output
  names(list.ts)=clim.var
  
  # compute the fun statistic (mean, min or max) for each point over the selected period
  if(!is.null(fun)) {
    
    FUN=get(fun)
    clim.fun=do.call(cbind,mclapply(as.list(clim.var),function(k) {apply(list.ts[[k]][,paste0('X',period)],1,FUN)},mc.cores=n.core))
    colnames(clim.fun)=clim.var
    clim.fun=data.frame(xy,clim.fun)
    
    if(!is.null(dir.out)) {
      
      nm.data=paste0(dir.out,"extraction_",min(period),'-',max(period),'_climate_',nrow(xy),"points",ifelse(!is.null(radius),paste0("_buffer",radius),""),"_eumedclim.RData")
      save(clim.fun,file=nm.data)
    }
    
    return(clim.fun)
    
  } else
    
    return(list.ts)
  
}


### Directory where climate data is going to be stored
dir.in='data/climate/'

### Output directory where extracted data are saved
dir.out='data/climate/'

### Data.frame of location coordinates
coord <- readRDS(file = "data/SiteProvCoord.RDS")

### Run 
ex=get.climate(xy=coord,clim.var="monthly",dir.in=dir.in,download=T,dir.out=dir.out,n.core=3)

# display output
ex

