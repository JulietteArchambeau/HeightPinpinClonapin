########################################################################################################################"
#                                                                                                                      #
#                    Checking distribution of heights in Caceres and Madrid                                            #
#                                       Figure S21                                                                     #
#                                                                                                                      #
#                                   Juliette Archambeau                                                                #
#                                       18/03/2022                                                                     #
#                                                                                                                      #
########################################################################################################################"

library(tidyverse) # CRAN v1.3.0
library(ggpubr)    # CRAN v0.2.1


# Function used in the script:
source("scripts/Functions/vir_lite.R") # available here: https://github.com/JulietteArchambeau/HeightPinpinClonapin/blob/master/scripts/Functions/vir_lite.R
ds=0.7 # parameter of the function vir_lite


data <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv") 

p <- ggarrange(
  data %>% dplyr::filter(site=="caceres") %>% 
    ggplot() + geom_histogram(aes(x=height),color="grey80", fill=vir_lite("pink",ds=ds),bins=50) + theme_bw() +
    labs(x="Height (mm)",y="Count"),
  
  data %>% dplyr::filter(site=="madrid") %>% 
    ggplot() + geom_histogram(aes(x=height),color="grey80", fill=vir_lite("deeppink",ds=ds),bins=50) + theme_bw() +
    labs(x="Height (mm)",y=""),
  
  nrow=1,labels=c("A) Caceres","B) Madrid"),font.label = list(size=20))


ggsave(p, file="figs/SuppInfo/ExploringData/DistributionHeightCaceresMadrid.png",width=14)
