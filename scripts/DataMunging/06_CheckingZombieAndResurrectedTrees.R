###############################################################################################"
##################                                                      #######################"
##################           Script from Juliette Archambeau            #######################"
##################                    July  2019                        #######################"
##################                                                      #######################"
##################               Function to check if I well            #######################"
##################            removed the resurrected and zombie        #######################"
##################                        trees                         #######################"
##################                                                      #######################"
###############################################################################################"


library(stringr)
library(tidyr)

data <- readRDS(file="data/AllDataPhenoClimSoil.RDS")




# A) Checking the resurrected trees    #####
############################################"


#### see function of the script "02_selecting_115_resurrected_trees_Portugal_Asturias.R"


sum.table <- data.frame(row.names = c("0","1","sum"))
res.trees <- data.frame()

for (i in c("bordeaux","asturias","portugal")){
  df <- data[data$site==i,]
  tab <-table(df$survival,df$age)
  tab <- rbind(tab,apply(tab,2,sum))
  tab <- as.data.frame(tab)
  rownames(tab) <- c("0","1","sum")
  colnames(tab) <- paste0(colnames(tab),sep="_",i)
  sum.table <- cbind(sum.table,tab)
  
  
  df$num.obs <- str_sub(df$obs, -1)
  df <- df[,c("site","block","clon","tree","num.obs","survival")]
  df <- spread(df,key="num.obs",value="survival")
  
  if(i=="portugal"){
    
    df <- df[!(df$`1`==1&df$`2`==1&df$`3`==1&df$`4`==0),]
    df <- df[!(df$`1`==1&df$`2`==1&df$`3`==1&df$`4`==1),]
    df <- df[!(df$`1`==0&df$`2`==0&df$`3`==0&df$`4`==0),]
    df <- df[!(df$`1`==1&df$`2`==1&df$`3`==0&df$`4`==0),]
    df <- df[!(df$`1`==1&df$`2`==0&df$`3`==0&df$`4`==0),]
    res.trees$`4` <- NA
    res.trees <- rbind(res.trees,df)
  } else{
    df <- df[!(df$`1`==1&df$`2`==1&df$`3`==1),]
    df <- df[!(df$`1`==0&df$`2`==0&df$`3`==0),]
    df <- df[!(df$`1`==1&df$`2`==1&df$`3`==0),]
    df <- df[!(df$`1`==1&df$`2`==0&df$`3`==0),]
    res.trees <- rbind(res.trees,df)
    
  }
  
}
res.trees <- res.trees[complete.cases(res.trees$site),]
sum.table
res.trees


## It's ok !! :) 



# B) Checking the zombie trees    #####
#######################################"

sum.table <- data.frame(row.names = c("0","1","sum"))
res.trees <- data.frame()

for (i in c("bordeaux","asturias","portugal")){
  df <- data[data$site==i,]
  tab <-table(df$survival,df$age)
  tab <- rbind(tab,apply(tab,2,sum))
  tab <- as.data.frame(tab)
  rownames(tab) <- c("0","1","sum")
  colnames(tab) <- paste0(colnames(tab),sep="_",i)
  sum.table <- cbind(sum.table,tab)
  
  
  df$num.obs <- str_sub(df$obs, -1)
  df <- df[,c("site","block","clon","tree","num.obs","survival")]
  df <- spread(df,key="num.obs",value="survival")
  
  if(i=="portugal"){
    
    df1 <- df[df$`1`==0&df$`2`==0&df$`3`==0&df$`4`==0,]
    df2 <- df[df$`1`==1&df$`2`==1&df$`3`==0&df$`4`==0,]
    df3 <- df[df$`1`==1&df$`2`==0&df$`3`==0&df$`4`==0,]
    res.trees$`4` <- NA
    res.trees <- rbind(res.trees,df1,df2,df3)
  } else{
    df1 <- df[df$`1`==0&df$`2`==0&df$`3`==0,]
    df2 <- df[df$`1`==1&df$`2`==0&df$`3`==0,]
    res.trees <- rbind(res.trees,df1,df2)
    
  }
  
}
res.trees <- res.trees[complete.cases(res.trees$site),]
sum.table
res.trees
unique(res.trees$site)

## It's ok !! :) 
