########################################################################################################################"
#                                                                                                                      #
#                            Correlations among climatic variables                                                     #
#                                                                                                                      #
#                                   Juliette Archambeau                                                                #
#                                       18/03/2022                                                                     #
#                                                                                                                      #
########################################################################################################################"

library(tidyverse)  # CRAN v1.3.0
library(ggcorrplot) # CRAN v0.1.3
library(ggbiplot)   # [github::vqv/ggbiplot] v0.55

data <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv")


# In the provenances ####


d.prov.selected <- data %>% 
  dplyr::select(starts_with("Q"),bio1_prov,bio5_prov,bio12_prov,bio14_prov) %>% 
  drop_na() %>% 
  distinct() %>% 
  scale() %>% 
  as_tibble() %>% 
  rename_at(vars(ends_with("_prov")), funs(str_sub(., 1, -6))) %>% 
  dplyr::rename("max.temp"=bio5,"mean.temp"=bio1,"mean.pre"=bio12,"min.pre"=bio14,
                "Northern Africa"=Q1, "Corsica"=Q2, "Central-Spain"=Q3, "French Atlantic"=Q4, "Iberian Atlantic"=Q5, "South-eastern Spain"=Q6)


# >>> Figure S7. Correlation plot ####

corr <- cor(d.prov.selected)

p <- ggcorrplot(corr, hc.order = TRUE, type = "lower",
                outline.col = "white",lab=T,
                legend.title = "Correlation\ncoefficients",
                tl.cex = 16) + 
  theme(legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        legend.position = c(0.2,0.8))

ggsave(p,file="figs/SuppInfo/ExploringData/CorrplotProv.png",height=10,width=12)


# >>> Figure S8. PCA ####

pca.prov.selected <- prcomp(d.prov.selected)

p <- ggbiplot(pca.prov.selected,varname.size =5,alpha = 0.3,labels.size = 10,
              varname.adjust = 1.5) +  ylim(-2.5, 1.6) +    xlim(-2.3, 2) + 
  theme_bw() +
  theme(plot.title = element_text(size=18),
        axis.title = element_text(size=18),
        axis.text = element_text(size=12))

ggsave(p,file="figs/SuppInfo/ExploringData/PcaProv.png",height=10,width=12)



# In the test sites ####

d.site <- data %>% 
  dplyr::select("age",
                "pre_summer_min_site",
                "pre_mean_1y_site",
                "tmn_min_1y_site",
                "tmx_max_1y_site",
                "pre_max_1y_site",
                "tmx_mean_1y_site") %>%
  drop_na() %>% 
  distinct()%>% 
  scale() %>% 
  as_tibble() %>% 
  rename_at(vars(ends_with("1y_site")), funs(str_sub(., 1, -9))) %>% 
  rename_at(vars(ends_with("_site")), funs(str_sub(., 1, -6))) %>% 
  dplyr::rename(mean.tmax=tmx_mean, 
                max.tmax=tmx_max, 
                min.presummer=pre_summer_min,
                min.tmn=tmn_min,
                mean.pre=pre_mean,
                max.pre=pre_max)


# >>> Figure S4. Correlation plot ####

cor <- cor(d.site)

p <- ggcorrplot(cor, 
                hc.order = TRUE, 
                type = "lower",
                outline.col = "white",
                lab=T,
                legend.title = "Correlation\ncoefficients",
                tl.cex = 16
) + 
  theme(legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        legend.position = c(0.2,0.8))

ggsave(p,file="figs/SuppInfo/ExploringData/CorrplotSite.png",height=10,width=12)


# >>> Figure S5. PCA ####

pca.site <- prcomp(d.site)

p <- ggbiplot(pca.site,varname.size =5,labels.size = 10,
              varname.adjust = 1.5) +   ylim(-2, 2) +    xlim(-3, 3) + 
  theme_bw() + geom_point(size=3,alpha = 0.5) +
  theme(plot.title = element_text(size=18),
        axis.title = element_text(size=18),
        axis.text = element_text(size=12))

ggsave(p,file="figs/SuppInfo/ExploringData/PcaSite.png",height=10,width=12)
