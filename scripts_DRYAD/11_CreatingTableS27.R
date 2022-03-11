########################################################################################################################"
#                                                                                                                      #
#                            Table S27 of the Supp Info in the Am Nat manuscript                                       #
#                                                                                                                      #
#                                   Juliette Archambeau                                                                #
#                                       11/03/2022                                                                     #
#                                                                                                                      #
########################################################################################################################"

library(brms)
library(dplyr)

mod <- readRDS(file="outputs/models/P1/MOD5.rds")

clones <- paste0("clon",1:6)


for(i in clones[1:length(clones)]){
  for(j in clones[1:length(clones)]){
    if(i==j){next}
    hyp <- hypothesis(mod,paste0(i,"__Intercept<",j,"__Intercept"),class="sd")
    if(i=="clon1"&j=="clon2"){
      df <- as.data.frame(hyp$hypothesis)
    } else {
      df <- bind_rows(df,hyp$hypothesis)
    }
  }}


df


print(xtable(df, type = "latex",digits=3), file = paste0("tables/Posteriors/M5_CompareTotGenVar.tex"), include.rownames=FALSE,sanitize.text.function = function(x) {x})
