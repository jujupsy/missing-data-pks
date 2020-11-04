#setwd("C:/Users/Juli/Dropbox/Psychologie B.Sc/8. SS 20/BA/parameter_estimation")


# packages
library(pks)
library(Rcpp)
library(RcppEigen)
library(parallel)
library(pbapply)

# compile cpp Funktions
#sourceCpp("./funktionen_cpp/EMmissblim.cpp"); sourceCpp("./funktionen_cpp/EMimblim.cpp")
#source("./funktionen_cpp/BLIM_combi.R")



fit_models <- function(data, workers = 10){
# init cluster

if(Sys.info()[['sysname']] == "Windows"){
    cores <- (parallel::detectCores() / 2) - 2 # hyperthreading:  cores/2 = physical cores
    cat("Running on Windows system. Using PSOCK as parallel system.\n") 
    cl <- makeCluster(min(cores, workers), type = "PSOCK", 
                      outfile = paste(gsub(":", "-", Sys.time()),"_worker.log",sep=""))  # no forking for windows :(
    cat("Compiling C++ functions... ")
    clusterEvalQ(cl, {
        library(pks)
        library(Rcpp)
        library(RcppEigen)
        sourceCpp("./funktionen_cpp/EMmissblim.cpp")
        sourceCpp("./funktionen_cpp/EMimblim.cpp")
        sourceCpp("./funktionen_cpp/EMblim.cpp")                
        source("./funktionen_cpp/BLIM_combi.R")
    })
    cat("DONE \n")
} else {
    cores <- parallel::detectCores()/2 
    library(pks)
    library(Rcpp)
    library(RcppEigen)
    cat("Compiling C++ functions... ")
    sourceCpp("./funktionen_cpp/EMmissblim.cpp")
    sourceCpp("./funktionen_cpp/EMimblim.cpp")
    sourceCpp("./funktionen_cpp/EMblim.cpp")                
    source("./funktionen_cpp/BLIM_combi.R")
    cat("DONE \n")
    cat("Running on Unix-like system. Using FORCK as parallel system.\n") 
    cl <- makeCluster(min(cores, workers), type = "FORK",
                      outfile = paste(gsub(":", "-", Sys.time()),"_worker.log",sep=""))
}
     

#load data
#WK <- "tKS"
#NinSet  <- 1000
#NofSets <- 1000
#fnames <- l(paste0("sim_data/data_mcar_", WK, "_", NinSet, "_", NofSets, ".rda"),
#             paste0("sim_data/data_ks_mnar_", WK, "_", NinSet, "_", NofSets, ".rda"),
#             paste0("sim_data/data_iks_mnar_", WK, "_", NinSet, "_", NofSets, ".rda"))
#f_names <- data$fnames
data_raw <- data
seed <- data$seed
files <- dir(data$fdir, pattern = ".rda", full.names = TRUE)

(t_start <- Sys.time())
data_estim <- NULL

for(mdl in c( "MissBLIM", "BLIM_COMP", "BLIM_MAW", "IMBLIM")){
    for(file in files){

       cat("Load next dataset.\n") 
       data <- readRDS(file)
       nMAX <- length(data)
       #data <- data[1] # test
       
       #next job is:
       cat(paste(mdl, ":", 
                   substr(data[[1]]$sim_condition, 1, 
                          nchar(data[[1]]$sim_condition) - 3), "\n"))
       
       if(Sys.info()[['sysname']] == "Windows"){
            clusterExport(cl, c("data", "mdl"), envir = environment())
            rm(data) # remove after transferred to cluster 
            gc()     # garbage collection. BAD HABIT but RAM!!!
       }

       tf <- pblapply(cl = cl, X = 1:nMAX, FUN = function(i,...) {
            
            m <- BLIM(data[[i]]$tModel$K,
                      data[[i]]$R, 
                      mdl, 
                      maxiter = 10000, 
                      fdb = FALSE, 
                      idx_k = TRUE)

            if(mdl != "MissBLIM"){
                m$mu0 <- m$mu1 <- rep(NA, times = m$nitems)
            }
            data.frame(cond = data[[i]]$sim_condition,
                       model = m$model,
                       para_q = c(paste0("beta_", 1:m$nitems),
                                paste0("eta_", 1:m$nitems),
                                paste0("mu0_", 1:m$nitems),
                                paste0("mu1_", 1:m$nitems),
                                paste0("PK_", 1:m$nstates)),
                       para =  c(rep("beta", m$nitems),
                                 rep("eta", m$nitems),
                                 rep("mu0", m$nitems),
                                 rep("mu1", m$nitems),
                                 rep("PK", m$nstates)),
                       estim_value = c(m$beta, 
                                       m$eta, 
                                       m$mu0, 
                                       m$mu1, 
                                       m$P.K),
                       true_value = c(data[[i]]$tModel$beta,
                                      data[[i]]$tModel$eta,
                                      data[[i]]$tModel$mu0,
                                      data[[i]]$tModel$mu1,
                                      data[[i]]$tModel$P.K),
                       mean_dist = mean(m$sym_dist),
                       set_num = data[[i]]$set_number,
                       converged = m$converged,
                       iterations = m$iterations,
                       N = m$N
                       )
        
                    
        })
        tf <- do.call(rbind, tf)
        data_estim <- rbind(data_estim, tf) 
        if(Sys.info()[['sysname']] == "Windows"){
            cat("Remove old datasets.\n")
            clusterEvalQ(cl, {rm(data); gc()})
        }else{
            cat("Remove old datasets.\n")
            rm(data)
            gc()
        }
        
    }
}

stopCluster(cl)
cat("Cluster stopped.\n")
(t_stop <- Sys.time())

(duration <- t_stop - t_start)
print(duration)
print(table(data_estim$converged))
cat("\n\n")
fn <- paste0("fitted_models_", data_raw$WK, "_", data_raw$NinSet, "_", data_raw$NofSets,"_", data_raw$seed, ".rda")
save(data_estim, file = fn)
cat(paste("Parameters of fitted models saved to:", fn, "\n\n"))
}
