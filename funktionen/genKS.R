# Funktion generiert WS nach DeChiusole et al. (2015)

### Datensimulation basierend auf:
# |Q| = 25
# |K| = 500
# beta_q ~ unif(0, .1)
# eta_q  ~ unif(0, .1)
# P(K)   ~ unif(0, 1) -> nomralized

genKS <- function(kQ = 25, kK = 500, seed = 42){
    cat(paste0("Generating new KS with seed: ", seed, "\n")) 
    set.seed(seed)
    # Potenzmenge von Q:
    potQ <- expand.grid(rep(list(0:1), kQ), KEEP.OUT.ATTRS = F)
    # 2^ == nrow(potQ)                  # check size


    # sample Elements of potQd that are part of K, without {} and Q  
    idx <- sort(sample(2:(nrow(potQ) - 1), size = kK - 2, replace = FALSE)) 

    K <- potQ[c(1, idx, nrow(potQ)), ]
    rownames(K) <- 1:kK

    # Parameter beta and eta
    beta <- runif(kQ, 0, .1)
    eta  <- runif(kQ, 0, .1)

    # P(K)
    P.K <- runif(kK, 0, 1)
    P.K <- P.K / sum(P.K)
    sum(P.K)

    true_para <- list(K = as.matrix(K), P.K = P.K, beta = beta, eta = eta, nitems = kQ, nstates = kK)
    return(true_para)
}




