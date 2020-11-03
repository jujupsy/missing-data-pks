# Name: sim_missing.R
# Description: Funcion to simulate missing data in KST
#              Possible dependencies between missings and M:
#                   - MCAR: P(Mq = 1 | q in K) = P(Mq = 1 | q not in K) = Pmcar 
#                   - MNAR: P(Mq = 1 | q in K) = MUq; 
#                           P(Mq = 1 | q not in K) = MUq_
#
# Input:  --
# Output: --
#
# Created:  30/Jun/2020, Julian Mollenhauer
#
##########################################
# working directory
#setwd("C:/Users/Juli/Dropbox/Psychologie B.Sc/8. SS 20/BA/Parameterschätzung/MissBLIM")
# packages


simMissingsBlim <- function(N, true_model, add_tm = TRUE){
    K    <- true_model$K
    P.K  <- true_model$P.K
    beta <- true_model$beta
    eta  <- true_model$eta
    mu0  <- true_model$mu0
    mu1  <- true_model$mu1

    nitems <- ncol(K)
    nstates <- nrow(K)
    
    # check input
    msg <- c("Number of", "parameters is incorrect or columns of K do not represent items.")
    if(length(beta) != nitems)
        stop(paste(msg[1], "beta", msg[2]))
    if(length(eta) != nitems)
        stop(paste(msg[1], "eta", msg[2]))
    if(length(mu0) != nitems)
        stop(paste(msg[1], "mu0", msg[2]))
    if(length(mu1) != nitems)
        stop(paste(msg[1], "mu1", msg[2]))       
    if(length(P.K) != nstates)
        stop(paste(msg[1], "P.K parameters is incorrect.")) 


    ## R* = RS complete response pattern ~ BlIM

    # sample N times k in K with prob. P.K
    idx_k <- sample(1:nrow(K), size = N, prob = P.K, replace = TRUE)

    # prob. to give correct answer    
    p_q_in_R <- t(t(K[idx_k, ]) * (1 - beta) + t(1 - K[idx_k, ]) * eta)

    # apply prob. and get complete response pattern
    RS <- t(apply(p_q_in_R, 1, function(r) rbinom(length(r), 1, t(r))))
    
    R <- RS
    # decimation mechanism leads to missings
    for(i in 1:nrow(R)){
        k <- as.numeric(K[idx_k[i], ])
        p_miss <- k * mu1 + (1 - k) * mu0
        idx_na <- rbinom(ncol(R), 1, p_miss)
        R[i, c(which(idx_na == 1))] <- NA
    }

    if(add_tm){
        ll <- list(R = R, idx_k = idx_k, RS = RS, tModel = true_model)
    } else {
        ll <- list(R = R, idx_k = idx_k, RS = RS)
    }

    return(ll)
}
