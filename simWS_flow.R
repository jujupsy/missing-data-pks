# Name: simWS_flow.R
# Description: Simulationsablauf fuer kuenstliche WS mit N = 100000
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
save(tKS, file = paste0("tKS_", s, ".rda"))
data <- data_sim(tks = paste0("tKS_", s), NinSet = 1000 * 100, NofSets = 200, cut = 50)
fit_models(data, workers = 10)
print("DONE :)")
