library(brms)
library(dplyr)

mod <- readRDS(file="outputs/models/P1/MOD5.rds")

clones <- paste0("clon",1:6)
as.data.frame(hypothesis(mod,paste0(i,"__Intercept<",j,"__Intercept"),class="sd"))


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

# Doesn't work
# 
# df$Hypothesis <- str_replace(df$Hypothesis,"(clon1__Intercept)","$\\sigma_{A_{NA}}^{2}$")
# df$Hypothesis <- str_replace(df$Hypothesis,"(clon2__Intercept)","$\\sigma_{A_{C}}^{2}$")
# df$Hypothesis <- str_replace(df$Hypothesis,"(clon3__Intercept)","$\\sigma_{A_{CS}}^{2}$")
# df$Hypothesis <- str_replace(df$Hypothesis,"(clon4__Intercept)","$\\sigma_{A_{FA}}^{2}$")
# df$Hypothesis <- str_replace(df$Hypothesis,"(clon5__Intercept)","$\\sigma_{A_{IA}}^{2}$")
# df$Hypothesis <- str_replace(df$Hypothesis,"(clon6__Intercept)","$\\sigma_{A_{SES}}^{2}$")
# df

print(xtable(df, type = "latex",digits=3), file = paste0("tables/Posteriors/M5_CompareTotGenVar.tex"), include.rownames=FALSE,sanitize.text.function = function(x) {x})
