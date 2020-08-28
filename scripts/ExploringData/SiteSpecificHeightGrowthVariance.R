
# In the manuscript,
# In the discussion, section "Predictiong new observations and provenances"

# What is the height growth variance in each site for each date of measurement? 
# ============================================================================

data <- readRDS(file="data/AllDataPhenoClimSoil.RDS")

# removing NAs for population structure
data <- data[!(is.na(data$Q1)),]

# removing NAs for height
data <- data[!(is.na(data$height)),]



for(s in unique(data$site)){
  
  for(age in unique(data$age[data$site==s])){
    
    cat(s," - ",age, ": ",var(log(data$height[data$site==s&data$age==age])),"\n")
    
    
  }
    
    
}
