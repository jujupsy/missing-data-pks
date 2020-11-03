# Name: MissBlim.R
# Description: Estimates parameters of BLIM / IMBLIM and MissBLIM Model using EM-Alogrithm (De Chiusole et al. 2015)
#              Aufbau grob orientiert an blim() Funktion (library(pks))
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
library(pks) # as.pattern und as.binmat werden verwendet


BLIM <- function(K, R, model = c("BLIM_MAW", "IMBLIM", "MissBLIM", "BLIM_COMP"), tol = 1e-07 , maxiter = 1000, fdb = FALSE, idx_k = FALSE){
    # Function to estimate parameter of missBLIM Model using EM algorithm (De Chiusole et al. 2015)
    
    model <- match.arg(model)
    print(model)
    # check arguments
    K <- as.matrix(K)
    R <- as.matrix(R)

    if(ncol(K) != ncol(R)) stop("Matrix K and R must have the same number of columns")

    R_raw <- R
    if(model == "BLIM_COMP"){
        # Antwortmuster mit NA entfernen
        R <- na.omit(R)
        
        if(nrow(R) == 0){
            cat("All response pattern contain missings. BLIM_COMP not possible!\n")
            para <- list()
            para$nitems <- ncol(K) 
            para$nstates <- nrow(K) 
            para$npat <- 0
            para$N <- 0
            para$N.RM <- 0
            para$converged <- NA
            para$K_asses_idx <- NA
            para$sym_dist <- NA
            para$model <- model
            para$beta <- rep(NA,  para$nitems)
            para$eta <- rep(NA,  para$nitems)
            para$iterations <- NA
            para$maxiter <- NA
            para$PK.R <- NA
            para$P.K <- rep(NA,  para$nstates)
            return(para)
        }
        
        if(nrow(R) < 10) 
            cat(paste("Only", nrow(R),  "response pattern contain no missings!\n"))
 
    }else if(model == "BLIM_MAW"){
        # MAW Transformation
        R[is.na(R)] <- 0
    } else {
         # NA als 9 kodieren
        R[is.na(R)] <- 9 # dann funktionieren as.pattern() und as.binmat() ohne Modifikation 
    }
    
    R_all <- R
    N.RM  <- as.pattern(R,   freq = TRUE)
    R     <- as.binmat(N.RM, uniq = TRUE)

    N       <- sum(N.RM) # Stichprobengroesse
    nitems  <- ncol(K)   # Anzahl an Problemen/Aufgaben q in Q
    nstates <- nrow(K)   # Anzahl an Wissenszustaenden k in K
    npat <- nrow(R)      # Anzahl Antwortmuster == N

    ## Initiale Werte 
    P.K  <- rep(1/nstates, nstates) # initiale Wkt. fuer P(K)
    beta <- rep(0.1, nitems)        # initiale Werte fuer beta, eta und mu
    eta  <- rep(0.1, nitems)
        
    mu0  <- rep(0.1, nitems)        # mu_q_
    mu1  <- rep(0.1, nitems)        # mu_q
    
    M <- (R == 9) # M: missing pattern
    W <- (R == 0) # W: wrong pattern
    R <- (R == 1) # R: right pattern

    iter <- 1
    maxdiff <- 2 * tol
    
    # convert to double for c++ function
    mode(R) <- mode(W) <- mode(M) <- mode(K) <- mode(N.RM) <- mode(P.K) <- mode(beta) <- mode(eta) <- mode(mu0) <- mode(mu1) <- "double"
    

    #### C++ Funktion waehlen
    if(model == "IMBLIM"){

        # IMBLIM
        P.M <- matrix(1, npat, nstates) * 
              (as.pattern(M * 1, freq = TRUE)[as.pattern(M * 1)] / N)
        mode(P.M) <- "double"  
        para <- emIMBLIMcpp(R, W, M, K, N.RM, P.K, beta, eta, P.M, 
                            maxiter, tol, fdb)

     }else if (model == "MissBLIM"){
         
        # MissBLIM
        para <- emMissBLIMcpp(R, W, M, K, N.RM, P.K, eta, beta, mu0, mu1, 
                              maxiter, tol, fdb)
        
    }else{
        
        # BLIM
        para <- emBLIMcpp(R, W, K, N.RM, P.K, beta, eta, maxiter, tol, fdb)

    }

    
    para$model <- model

    # warning if maxiter was reached 
    if(para$maxiter) warning(paste("Maximum # of", maxiter, "was reached!"))
    if(!para$converged) warning(paste("EM-Algorithm did NOT converge"))
    
    # Assesement of KS
    K_asses <- apply(para$PK.R, 2, function(x) which(x == max(x))[1])
    para$K_asses_idx <- apply(R_all, 1, function(x) K_asses[which(names(N.RM) == paste(x, collapse = ""))])
    
    # symmetric Distance between K_asses and K if rownames of R correspond to idx of latent KS in K
    if(idx_k && model != "BLIM_COMP"){
        para$sym_dist <- apply(cbind(as.numeric(rownames(R_raw)), para$K_asses_idx), 1, 
                               function(x) sum(K[x[1], ] !=  K[x[2], ]))
    }else if(idx_k && model == "BLIM_COMP"){
        para$sym_dist <- apply(cbind(as.numeric(rownames(R_raw)[-attributes(R_all)$na.action]), para$K_asses_idx), 1, 
                               function(x) {sum(K[x[1], ] !=  K[x[2], ])})
        
    }

    # add information to output
    para$nitems <- nitems
    para$nstates <- nstates
    para$npat <- npat
    para$N <- N
    para$N.RM <- N.RM
    return(para)
}



