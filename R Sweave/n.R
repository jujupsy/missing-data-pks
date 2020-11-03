# Name: n.R
# Description: Angepasste Modelle für die verschiedenen Stichprobengroessen 
#              vorverarbeiten, da Rechendauer (einlesen etc.) zu lange dauert
#              und sonst das kompilieren der .Rnw Datei zu lange braucht
#
# Input:  fitted Models in data/
# Output: "data/fitted_models_n.rda" und "data/ag_bias_n.rda"
#
# Created:  02/Sep/2020, Julian Mollenhauer
#
##########################################
# working directory
setwd("C:/Users/Juli/Dropbox/Psychologie B.Sc/8. SS 20/BA/Schriftfassung")
# packages

if("fitted_models_n.rda" %in% dir("data/")){
    load("data/fitted_models_n.rda")
} else {
    data <- NULL
    ## bias related to number of subjects
    for(f in dir("data/n/", full.names = TRUE, pattern = "fitted_models")){
        load(f)
        data_estim$n <- max(data_estim$N)
        data <- rbind(data, data_estim)
    }
    save(data, file = "data/fitted_models_n.rda")
}


# Table Bias 
data$bias <- data$estim_value - data$true_value

data1 <- subset(data, para %in% c("beta", "eta"))

data1$model <- factor(data1$model)
data1$para  <- factor(data1$para)
# works also with subset of data!
level_cond <- c(paste0("mc_", (1:4)*10), paste0("ks_", (1:4) * 10), paste0("iks_C", 1:5))
data1$cond  <- factor(data1$cond, levels = level_cond[which(level_cond %in% levels(factor(data1$cond)))])

# only 10 and 20 percent missing und c2, c4
dat <- droplevels(data1[data1$cond %in% c("mc_10", "mc_20", "ks_10", "ks_20", "iks_C2", "iks_C4"), ])

# sample size 
n_comp <- aggregate(N ~ cond + model + n, dat[dat$model == "BLIM_COMP", ], mean)

# remove estimates that did not converge
#data_estim_raw <- data_estim

not_converged <- unique(dat[(!is.na(dat$converged) & dat$converged == FALSE), c("cond", "model", "N", "n", "iterations", "set_num")])
ftable(not_converged$model, not_converged$cond, not_converged$n)


dat<- na.omit(dat[dat$converged, ])

ag_bias <- aggregate(bias ~ cond + model + para + n, dat, mean, drop = FALSE)
ag_bias$SD <- aggregate(bias ~ cond + model + para + n, dat, sd, drop = FALSE)[, 5]

save(ag_bias, file = "data/ag_bias_n.rda")




