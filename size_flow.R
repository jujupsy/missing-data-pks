# Name: size_flow.R
# Description: Simulationsablauf verschiedene Stichprobengroessen mit kuenstlicher WS
#
# Input:  --
# Output: --
#
# Created:  29/Jul/2020, Julian Mollenhauer
#
##########################################
# working directory
# setwd("C:/Users/Juli/Dropbox/Psychologie B.Sc/8. SS 20/BA/simulations")
# packages
source("data_sim.R")
source("fit_models.R")
source("funktionen/genKS.R")

s <- 1432
tKS <- genKS(seed = s)
save(tKS, file = paste0("n_tKS_", s, ".rda"))

ns <- sort(c(2^c(0:8) * 100)) # SP Groesse

for(n in ns){
    data <- data_sim(tks = paste0("n_tKS_", s), NinSet = n, NofSets = 200, cut = 10)
    fit_models(data, workers = 10)
}
print("DONE :)")
