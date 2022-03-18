########################################################################################################################"
#                                                                                                                      #
#                            Details on the experimental design                                                        #
#                                                                                                                      #
#                                   Juliette Archambeau                                                                #
#                                       18/03/2022                                                                     #
#                                                                                                                      #
########################################################################################################################"

library(tidyverse) # CRAN v1.3.0


# Table S1   ####
# ========   "

# >> Number of observations in the entire dataset and in each partition

# A. In the entire dataset

data <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv")

tab <- as.data.frame(table(data$site))
colnames(tab) <-c("Level","All") 
tab$Level <- str_to_title(tab$Level)

data$site_age <- paste0(data$site, " - ", data$age, " months old")
data$site_age <- str_to_sentence(data$site_age)
tab_add <- as.data.frame(table(data$site_age))
colnames(tab_add) <-c("Level","All") 

tab <- bind_rows(tab,tab_add)

tab <- bind_rows(data.frame(Level="All sites",All=dim(data)[[1]] ), tab)


# B. In the P1 partition 

# B.a. Train dataset

data <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv") %>%  dplyr::filter(P1=="train")

tab_add <- as.data.frame(table(data$site))
colnames(tab_add) <-c("Level","Train1") 
tab_add$Level <- str_to_title(tab_add$Level)

data$site_age <- paste0(data$site, " - ", data$age, " months old")
data$site_age <- str_to_sentence(data$site_age)
tab_add_2 <- as.data.frame(table(data$site_age))
colnames(tab_add_2) <-c("Level","Train1") 

tab_add <- bind_rows(tab_add,tab_add_2)

tab_add <- bind_rows(data.frame(Level="All sites",Train1=dim(data)[[1]] ), tab_add)
tab <- left_join(tab,tab_add)
tab$Train1  <- tab$Train1 %>% replace_na(0)


# B.b. Test dataset

data <- read_csv("data_DRYAD/TestP1prepared.csv")

tab_add <- as.data.frame(table(data$site))
colnames(tab_add) <-c("Level","Test1") 
tab_add$Level <- str_to_title(tab_add$Level)

data$site_age <- paste0(data$site, " - ", data$age, " months old")
data$site_age <- str_to_sentence(data$site_age)
tab_add_2 <- as.data.frame(table(data$site_age))
colnames(tab_add_2) <-c("Level","Test1") 

tab_add <- bind_rows(tab_add,tab_add_2)

tab_add <- bind_rows(data.frame(Level="All sites",Test1=dim(data)[[1]] ), tab_add)
tab <- left_join(tab,tab_add)
tab$Test1  <- tab$Test1 %>% replace_na(0)


# C. In the P2 partition 

# C.a. Train dataset

data <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv") %>%  dplyr::filter(P2=="train")

tab_add <- as.data.frame(table(data$site))
colnames(tab_add) <-c("Level","Train2") 
tab_add$Level <- str_to_title(tab_add$Level)

data$site_age <- paste0(data$site, " - ", data$age, " months old")
data$site_age <- str_to_sentence(data$site_age)
tab_add_2 <- as.data.frame(table(data$site_age))
colnames(tab_add_2) <-c("Level","Train2") 

tab_add <- bind_rows(tab_add,tab_add_2)

tab_add <- bind_rows(data.frame(Level="All sites",Train2=dim(data)[[1]] ), tab_add)
tab <- left_join(tab,tab_add)
  tab$Train2  <- tab$Train2 %>% replace_na(0)


# C.b. Test dataset

data <- read_csv("data_DRYAD/TestP2prepared.csv")

tab_add <- as.data.frame(table(data$site))
colnames(tab_add) <-c("Level","Test2") 
tab_add$Level <- str_to_title(tab_add$Level)

data$site_age <- paste0(data$site, " - ", data$age, " months old")
data$site_age <- str_to_sentence(data$site_age)
tab_add_2 <- as.data.frame(table(data$site_age))
colnames(tab_add_2) <-c("Level","Test2") 

tab_add <- bind_rows(tab_add,tab_add_2)

tab_add <- bind_rows(data.frame(Level="All sites",Test2=dim(data)[[1]] ), tab_add)
tab <- left_join(tab,tab_add)
tab$Test2  <- tab$Test2 %>% replace_na(0)


# D. In the P3 partition 

# D.a. Train dataset

data <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv") %>%  dplyr::filter(P3=="train")

tab_add <- as.data.frame(table(data$site))
colnames(tab_add) <-c("Level","Train3") 
tab_add$Level <- str_to_title(tab_add$Level)

data$site_age <- paste0(data$site, " - ", data$age, " months old")
data$site_age <- str_to_sentence(data$site_age)
tab_add_3 <- as.data.frame(table(data$site_age))
colnames(tab_add_3) <-c("Level","Train3") 

tab_add <- bind_rows(tab_add,tab_add_3)

tab_add <- bind_rows(data.frame(Level="All sites",Train3=dim(data)[[1]] ), tab_add)
tab <- left_join(tab,tab_add)
tab$Train3  <- tab$Train3 %>% replace_na(0)


# D.b. Test dataset

data <- read_csv("data_DRYAD/TestP3prepared.csv")

tab_add <- as.data.frame(table(data$site))
colnames(tab_add) <-c("Level","Test3") 
tab_add$Level <- str_to_title(tab_add$Level)

data$site_age <- paste0(data$site, " - ", data$age, " months old")
data$site_age <- str_to_sentence(data$site_age)
tab_add_3 <- as.data.frame(table(data$site_age))
colnames(tab_add_3) <-c("Level","Test3") 

tab_add <- bind_rows(tab_add,tab_add_3)

tab_add <- bind_rows(data.frame(Level="All sites",Test3=dim(data)[[1]] ), tab_add)
tab <- left_join(tab,tab_add)
tab$Test3  <- tab$Test3 %>% replace_na(0)



row.names(tab) <- tab$Level
tab$Level <- NULL
tab <- round(tab,0)

# Table S1:
print(xtable(tab, type = "latex",digits=0), file = "tables/ExperimentalDesign/ExpDesignNbObs.tex")



# Table S2   ####
# ========   "

# Provenance information: provenance codes used in the study, number of genotypes, trees and observations
# (an observation being a height-growth measurement in a given year on one individual) per provenance.

data <-  read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv")


dfclon <- data %>%  dplyr::select(prov,clon) %>% distinct() %>% group_by(prov) %>% dplyr::count()
colnames(dfclon) <- c("Provenances","Number of clones")
dftree <- data %>%  dplyr::select(prov,tree) %>% distinct() %>% group_by(prov) %>% dplyr::count()
colnames(dftree) <- c("Provenances","Number of trees")

dfobs <- data %>%  group_by(prov) %>% tally()
colnames(dfobs) <- c("Provenances","Number of observations")

df <- inner_join(dfclon,dftree)
df <- inner_join(df,dfobs)


print(xtable(df, type = "latex",digits=0), 
      file = "tables/ExperimentalDesign/ExpDesignProvenances.tex", 
      include.rownames=FALSE)


# Table S3   ####
# ========   "

# Mean proportion belonging to each gene pool for each provenance.

# Northern africa gene pool => NA => Q1
# Corsican gene pool => C => Q2
# Central Spain gene pool => CS => Q3
# French Atlantic gene pool => FA => Q4
# Iberian Atlantic gene pool => IA => Q5
# South-eastern Spain gene pool => SES => Q6

data <- read_csv("data_DRYAD/HeightClimateSoilData_33121obs_32variables.csv")

ps <- data %>%
  group_by(prov) %>% 
  summarise_at(vars(paste0(rep("Q",6),1:6)), mean)

colnames(ps) <- c("Provenance","NA","C","CS","FA","IA","SES")

print(xtable(ps, type = "latex",digits=3), 
      file = "tables/ExperimentalDesign/PopStructureByProv.tex", 
      include.rownames=FALSE)

