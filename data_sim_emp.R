# Name: data_sim_emp.R
# Description: 
#       Simulate data for pks following the procedure from 
#           De Chiusole et al. (2015) from tKS_fpi from FPI-R
#       mu_q und mu_q_ are adjusted (smaller compared to sim WS)
#
# Input:
#     - functionen/simMissingsBlim.R
#     - tKS.rda 
#        # true knowledge strucuter simulated according to 
#           De Chiusole et al. (2015)
#           
# Output: 
#     - data_mcar.rda
#     - data_ksmnar.rda
#     - data_iksmnar.rda
#
# Created:  06/Jul/2020, Julian Mollenhauer
#
##########################################
# working directory
#setwd("C:/Users/Juli/Dropbox/Psychologie B.Sc/8. SS 20/BA/parameter_estimation")

# packages and functions
library(pks)
source("./funktionen/simMissingsBlim.R")




data_sim_emp <- function(NinSet, NofSets, tks = "tKS", seed = 42, cut = 100){

# load true knowledge structure:
# - always same KS (tKS)
#   - |Q| = 25, |K| = 500
# - always same beta eta but mu0 and mu1 according to missing generating process
load_tks <- T
set.seed(seed)
cat(paste("Seed for data simulation is:", seed, "\n\n"))
if(load_tks){
    cat(paste0("load KS from file ", tks, ".rda\n"))
    load(paste0(tks, ".rda"))
    WK <- tks
    
} else {
    cat("use KS from pks\n")
    data(chess)
    WK <- "dst1"
    K <- chess$dst1
    N.R <- chess$N.R
    tKS <- blim(K, N.R, method = "MDML")
}



# sample size in each simulated dataset

N <- NinSet
n_datasets <- NofSets
#cut <- 100 

cat(paste("Sample size in each dataset is: ", N, "\n"))

# create new folder in sim_data/
path <- paste0("sim_data/", WK, "_", N, "_", n_datasets, "_s", seed, "_c", cut)
dir.create(path, 
           showWarnings = F)
if(length(dir(path, pattern = ".rda", full.names = TRUE)) > 3){
    cat(paste("Dataset with seed", seed, "already exists in", path, ":)\n\n"))
    return(list(fdir = path, WK = WK, NinSet = NinSet, NofSets = NofSets, seed = seed))
}else{
    cat(paste("Dataset with seed", seed, 
              "does not exist yet. Starting to simulate the data... :)\n"))
}
### Missing Completely at Random (MCAR) ###
# P(Mq = 1 | q in K) = P(Mq = 1 | q not in K) = Pmcar 

# number of datasets for each mu0 mu1 combination (each row in mcar_mu)
n_mcar <- n_datasets
pmcar <- c(.01, .05, .1, .15)
cond_mcar <- paste0("mc_", pmcar*100)
mcar_mu <- list(mu0 = pmcar, mu1 = pmcar)

data_mcar <- list()
idx = 1
for(i in seq_along(pmcar)){
    for(j in seq_len(n_mcar)){
        if((idx %% cut) == 0){
            fn1 <- paste0(path, "/data_mcar_", WK, "_", N, "_", 
                          n_datasets, "_", seed,"_#", idx, ".rda")
            saveRDS(data_mcar, file = fn1)
            data_mcar <- list()
        }
        tKS$mu0 <- rep(mcar_mu$mu0[i], tKS$nitems)
        tKS$mu1 <- rep(mcar_mu$mu1[i], tKS$nitems)
        t <- list(simMissingsBlim(N, true_model = tKS, add_tm = TRUE))
        t[[1]]$sim_condition <- cond_mcar[i]
        t[[1]]$set_number <- j
        names(t) <- paste0("mcar_", i, "_", j)
        data_mcar <- append(data_mcar, t)
        idx <- idx + 1
    }
}
fn1 <- paste0(path, "/data_mcar_", WK, "_", N, "_", 
              n_datasets, "_", seed,"_#", idx, ".rda")
saveRDS(data_mcar, file = fn1)
data_mcar <- list()

cat(paste("MCAR simulated and saved to", path, "\n"))

### ks-Missing Not at Random (ks-MNAR) ###
# ks-MNAR: P(Mq = 1 | q in K) = 0 
#          P(Mq = 1 | q not in K) = pks


# number of datasets for each mu0 mu1 combination (each row in mcar_mu)
n_ksmnar <- n_datasets
p_ks <- c(.02, .10, .20, .30)
cond_ks <- paste0("ks_", p_ks*50)
ksmnar_mu <- list(mu1 = rep(0, length(p_ks)), mu0 = p_ks)

data_ks_mnar <- list()
idx = 1
for(i in seq_along(p_ks)){
    for(j in seq_len(n_ksmnar)){
        if((idx %% cut) == 0){
                fn2 <- paste0(path, "/data_ks_mnar_", WK, "_", N, "_", 
                              n_datasets, "_", seed,"_#", idx, ".rda")
                saveRDS(data_ks_mnar, file = fn2)
                data_ks_mnar <- list()
        }
        tKS$mu0 <- rep(ksmnar_mu$mu0[i], tKS$nitems)
        tKS$mu1 <- rep(ksmnar_mu$mu1[i], tKS$nitems)
        t <- list(simMissingsBlim(N, true_model = tKS, add_tm = TRUE))
        t[[1]]$sim_condition <- cond_ks[i]
        t[[1]]$set_number <- j
        names(t) <- paste0("ksmnar_", i, "_", j)
        data_ks_mnar <- append(data_ks_mnar, t)
        idx <- idx + 1
    }
}
fn2 <- paste0(path, "/data_ks_mnar_", WK, "_", N, "_", 
              n_datasets, "_", seed,"_#", idx, ".rda")
saveRDS(data_ks_mnar, file = fn2)
data_ks_mnar <- list()

cat(paste("ks-MNAR simulated and saved to", path, "\n"))

### iks-Missing Not at Random (ks-MNAR) ###
# ks-MNAR: P(Mq = 1 | q in K) = mu1 
#          P(Mq = 1 | q not in K) = mu0


# number of datasets for each mu0 mu1 combination (each row in mcar_mu)
n_iksmnar <- n_datasets

# sample mu: mu ~ unif(a, b)

a_mu1 <- rep(c(.10, .05, .01, 0), each = n_iksmnar) 
b_mu1 <- rep(c(.15, .10, .05,.01), each = n_iksmnar) 
cond_iks <- rep(paste0("iks_C", 1:4), each = n_iksmnar)
set_number <- rep(1:n_iksmnar, times = 4)
a_mu0 <- rev(a_mu1)
b_mu0 <- rev(b_mu1)

data_iks_mnar <- list()
idx <- 1
for(i in seq_along(a_mu1)){
    if((idx %% cut) == 0){
        fn3 <- paste0(path, "/data_iks_mnar_", WK, "_", N, "_", 
                      n_datasets, "_", seed,"_#", idx, ".rda")
        saveRDS(data_iks_mnar, file = fn3)
        data_iks_mnar <- list()
    }
    tKS$mu0 <- runif(tKS$nitems, a_mu0[i], b_mu0[i])
    tKS$mu1 <- runif(tKS$nitems, a_mu1[i], b_mu1[i])
    t <- list(simMissingsBlim(N, true_model = tKS, add_tm = TRUE))
    t[[1]]$sim_condition <- cond_iks[i]
    t[[1]]$set_number <- set_number[i]
    names(t) <- paste0(cond_iks[i], "_", set_number[i])
    data_iks_mnar <- append(data_iks_mnar, t)
    idx <- idx + 1
}
fn3 <- paste0(path, "/data_iks_mnar_", WK, "_", N, "_", 
              n_datasets, "_", seed,"_#", idx, ".rda")
saveRDS(data_iks_mnar, file = fn3)

cat(paste("iks-MNAR simulated and saved to", path, "\n"))
cat("Data simulation is finished.\n\n")
return(list(fdir = path, WK = WK, NinSet = NinSet, NofSets = NofSets, seed = seed))
}
