# Correlation between PEAs and proportion of belonging to each gene pool

# 06/05/2020

mod <- readRDS(file="outputs/models/P1/MOD7.rds") 
data <- mod$data
sub <- data[,c("prop_Q1","prop_Q2","prop_Q3","prop_Q4","prop_Q5","prop_Q6","count_all.sc")]
sub <- unique(sub)
cor(sub)


mod <- readRDS(file="outputs/models/P1/MOD8.rds") 
data <- mod$data
sub <- data[,c("prop_Q1","prop_Q2","prop_Q3","prop_Q4","prop_Q5","prop_Q6","count_all_350.sc")]
sub <- unique(sub)
cor(sub)
