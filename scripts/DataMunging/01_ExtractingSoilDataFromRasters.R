############################################################################################"
##################  Script from Juliette Archambeau - November 2018  #######################"
##################   Setting a common db with coord and soil data    #######################"
############################################################################################"

# libraries
library(raster)

# dataframe of location coordinates
coord <- readRDS(file = "data/SiteProvCoord.RDS")

# dataset where we are going to add soil data
df <- coord

###########################################################################"
# prov CEN = 0 for all columns
# We are going to calculate the mean of the four pixels around prov CEN.
coord.CEN <- as.data.frame(matrix(NA, 4,2, dimnames = list(c("pixel.left","pixel.right","pixel.above","pixel.below"),
                                                        c("longitude","latitude"))))
coord.CEN["pixel.left",] <- c(-4.5050,40.27836)
coord.CEN["pixel.right",] <- c(-4.4813,40.27836)
coord.CEN["pixel.above",] <- c(-4.4912807,40.2893)
coord.CEN["pixel.below",] <- c(-4.4912807,40.2646)
##########################################################################"


########  1- depth roots
depth <- raster("data/soil/DERIVED_LAYERS/STU_EU_DEPTH_ROOTS_WGS84.tif")
df$depth_roots <- raster::extract(depth, coord)
df["CEN","depth_roots"] <- mean(raster::extract(depth, coord.CEN))
remove(depth)

########  2- clay content - topsoil
clay_top <- raster("data/soil/DERIVED_LAYERS/STU_EU_T_CLAY_WGS84.tif")
df$clay_top <- raster::extract(clay_top, coord)
df["CEN","clay_top"] <- mean(raster::extract(clay_top, coord.CEN))
remove(clay_top)

########  3- silt content - topsoil
silt_top <- raster("data/soil/DERIVED_LAYERS/STU_EU_T_SILT_WGS84.tif")
df$silt_top <- raster::extract(silt_top, coord)
df["CEN","silt_top"] <- mean(raster::extract(silt_top, coord.CEN))
remove(silt_top)

########  4- silt content - topsoil
sand_top <- raster("data/soil/DERIVED_LAYERS/STU_EU_T_SAND_WGS84.tif")
df$sand_top <- raster::extract(sand_top, coord)
df["CEN","sand_top"] <- mean(raster::extract(sand_top, coord.CEN))
remove(sand_top)

########  5- Total available water content from PTF - topsoil 
water_top <- raster("data/soil/DERIVED_LAYERS/STU_EU_T_TAWC_WGS84.tif")
df$water_top <- raster::extract(water_top, coord)
df["CEN","water_top"] <- mean(raster::extract(water_top, coord.CEN))
remove(water_top)

########  6- clay content - subsoil
clay_sub <- raster("data/soil/DERIVED_LAYERS/STU_EU_S_CLAY_WGS84.tif")
df$clay_sub <- raster::extract(clay_sub, coord)
df["CEN","clay_sub"] <- mean(raster::extract(clay_sub, coord.CEN))
remove(clay_sub)

########  7- silt content - subsoil
silt_sub <- raster("data/soil/DERIVED_LAYERS/STU_EU_S_SILT_WGS84.tif")
df$silt_sub <- raster::extract(silt_sub, coord)
df["CEN","silt_sub"] <- mean(raster::extract(silt_sub, coord.CEN))
remove(silt_sub)

########  8- silt content - subsoil
sand_sub <- raster("data/soil/DERIVED_LAYERS/STU_EU_S_SAND_WGS84.tif")
df$sand_sub <- raster::extract(sand_sub, coord)
df["CEN","sand_sub"] <- mean(raster::extract(sand_sub, coord.CEN))
remove(sand_sub)

########  9- Total available water content from PTF - subsoil 
water_sub <- raster("data/soil/DERIVED_LAYERS/STU_EU_S_TAWC_WGS84.tif")
df$water_sub <- raster::extract(water_sub, coord)
df["CEN","water_sub"] <- mean(raster::extract(water_sub, coord.CEN))
remove(water_sub)


####################################################"
# In the sub layers, there some O that are actually NAs. 
df[df$water_sub==0  &df$clay_sub==0 & df$silt_sub==0& df$sand_sub==0,
   c("clay_sub","sand_sub","silt_sub","water_sub")] <- NA
####################################################"

### Separate site and prov data
provs <- df[c(1:35),]
colnames(provs) <- paste(colnames(provs),"prov",sep="_") 
saveRDS(provs,file="data/ProvSoilData.RDS")

sites <- df[c(36:40),]
colnames(sites) <- paste(colnames(sites),"site",sep="_") 
saveRDS(sites,file="data/SiteSoilData.RDS")

####################################################################"
##############           COMMENTS         ##########################"

# 16 locations where sib layers = 0 => better to choose top layers. 
