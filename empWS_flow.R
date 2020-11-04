# Name: empWS_flow.R
# Description: Simulationsablauf empirische WS aus dem FPI-R hergeleitet
#
# Input:  --
# Output: --
#
# Created:  29/Aug/2020, Julian Mollenhauer
#
##########################################
# working directory
# setwd("C:/Users/Juli/Dropbox/Psychologie B.Sc/8. SS 20/BA/simulations")
# packages
source("data_sim_emp.R")
source("fit_models.R")

#########
N <- 3700
# FPI-R Gesundheitssorgen
data <- data_sim_emp(tks = "tKS_fpi_ges", NinSet = N, NofSets = 200, cut = 50)
fit_models(data, workers = 10)

# FPI-R Gehemmtheit
data <- data_sim_emp(tks = "tKS_fpi_geh", NinSet = N, NofSets = 200, cut = 50)
fit_models(data, workers = 10)

# FPI-R Erregbarkeit
data <- data_sim_emp(tks = "tKS_fpi_err", NinSet = N, NofSets = 200, cut = 50)
fit_models(data, workers = 10)

#########
N <- 10000
# FPI-R Gesundheitssorgen
data <- data_sim_emp(tks = "tKS_fpi_ges", NinSet = N, NofSets = 200, cut = 50)
fit_models(data, workers = 10)

# FPI-R Gehemmtheit
data <- data_sim_emp(tks = "tKS_fpi_geh", NinSet = N, NofSets = 200, cut = 50)
fit_models(data, workers = 10)

# FPI-R Erregbarkeit
data <- data_sim_emp(tks = "tKS_fpi_err", NinSet = N, NofSets = 200, cut = 50)
fit_models(data, workers = 10)

print("DONE :)")

