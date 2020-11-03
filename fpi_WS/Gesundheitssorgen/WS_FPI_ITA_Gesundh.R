# Name: WS_FPI_ITA.txt
# Description: generate Knowledgestructure from FPI-R scale "Gesundheitssorgen"
#
# Input:  --
# Output: --
#
# Created:  07/Sep/2020, Julian Mollenhauer
#
##########################################
# working directory
setwd("C:/Users/Juli/Dropbox/Psychologie B.Sc/8. SS 20/BA/simulations/fpi_WS/Gesundheitssorgen")

# packages
library("pks")
library("relations")
# BiocManager::install("Rgraphviz")
library("Rgraphviz")

# read data and code 999 as NA 
raw <- read.table("../fpi_data/fgjn99fr19_pd.txt", header = TRUE, na.strings = 999)

# missing data table
table(rowSums(is.na(raw[, 1:138])))

# Gesundheitssorgen: items and coding
item <- c(10, 18, 31, 38, 62, 65, 68, 70, 84, 89, 117, 127)
ge <- raw[, abs(item)]
names(ge)
# 2 -> 0 ("stimmt nicht"); 1 = "stimmt"
ge[ge == 2] <- 0

# flip  reversed poled items and remove patterns containing NA
if(any(item < 0))
    ge[, which(item < 0)] <- (1 -  ge[,  which(item < 0)]) 

nrow(ge) # 3740
table(rowSums(is.na(ge)))
#    0    1 
# 3706   34 
ge <- na.omit(ge)
nrow(ge) #  3706

# sample information
SP <- raw[-attributes(ge)$na.action, c("geschlecht", "alter")]
mean(SP$alter)  # 45.81
sd(SP$alter)    # 17.65
range(SP$alter) # 16 95

# 1: weiblich, 2:maennlich
table(SP$geschlecht)
# percentage female
mean(SP$geschlecht == 1) # 0.466

# distinct patterns
nrow(unique(ge)) # 1350

### ITA Method (Schrepp, 1999, Leeuve, 1974)

# b values, precedence relation violated 
b <- matrix(0, ncol(ge), ncol(ge), dimnames = list(names(ge), names(ge)))
for(i in 1:ncol(ge)){
    for(j in 1:ncol(ge)){
        b[j, i] <- sum(ge[, i] == 1 & ge[,j] == 0)
    }
}

# Schwelle: L = 277
 L <- 277
rel <- relation(incidence = b <= L)
relation_violations(rel, property = "transitive") # 0 :D

prece <- relation_incidence(rel)

K <- as.pattern(t(prece), as.set = TRUE)
K <- closure(K, operation = "intersection")
K <- closure(K, operation = "union")

K <- as.binmat(K)

nrow(K) # 560


# rename
idx <- order(colnames(prece))
prece <- prece[idx, idx]
lab.items <- paste("  I", sort(abs(item)), "  ")
dimnames(prece) <- list(lab.items, lab.items)
rel <- relation(incidence = prece)

# Surmise Relation plotten (in LATEX)
#save(rel, file = "precedence_fpi_ges.rda")
#load("precedence_fpi_ges.rda")
plot(rel, main = "",
     attrs = list(graph = list(rankdir = "BT"),
                  edge = list(arrowsize = NULL),
                  node = list(shape = "ellipse",
                  fixedsize = FALSE)))

## Properties of Rel
relation_violations(rel, property = "reflexiv")       # 0 -> reflexiv
relation_violations(rel, property = "transitive")     # 0 -> transitiv
relation_violations(rel, property = "antisymmetric")  # 0 -> antisymmetrisch --> well-graded (Halbordnung auf der Menge der Items)

#library("kst") # double-check properties
#kstr <- kstructure(as.pattern(K, as.set = T))
#kstructure_is_wellgraded(kstr) # ->



## parameter estimation 
# load if saved estimates are in directory
if("estimates_fpi_ges.rda" %in% dir()){
    load("estimates_fpi_ges.rda")
} else{
    # Minimum Discrepency
    blim.md <- blim(K,  table(as.pattern(ge)), method = "MDML")

    # Maximum Likelihood
    library(pks)
    library(Rcpp)
    library(RcppEigen)
    #sourceCpp("../funktionen_cpp/EMmissblim.cpp")
    #sourceCpp("../funktionen_cpp/EMimblim.cpp")
    sourceCpp("../../funktionen_cpp/EMblim.cpp")                
    source("../../funktionen_cpp/BLIM_combi.R")

    blim.ml  <- BLIM(K, ge, "BLIM_COMP", maxiter = 10000, fdb = TRUE, idx_k = FALSE)
    save(blim.md, blim.ml, file = "estimates_fpi_ges.rda")

}


#compare md to ml estimation
est_beta <- as.matrix(cbind(blim.md$beta, blim.ml$beta))

est_eta <- as.matrix(cbind(blim.md$eta, blim.ml$eta))

colnames(est_beta) <- colnames(est_eta) <- c("md", "ml")

par(mfrow = c(1, 2))
barplot(t(est_beta), beside = T, legend.text = T, main = "beta")
barplot(t(est_eta),  beside = T, legend.text = T, main = "eta")


# save ks for simulation:
#rownames K
K <- as.matrix(K)
rownames(K) <- 1:nrow(K)
tKS <- list(K = K, P.K = blim.ml$P.K, beta = blim.ml$beta, eta = blim.ml$eta, nitems = ncol(ge), nstates = nrow(K))
save(tKS, file = "tKS_fpi_ges.rda")


# Identifizierbarkeit
blim.ml$N.R <- blim.ml$N.RM
blim.ml$K <- K
jac <- jacobian(blim.ml, blim.ml$P.K, blim.ml$beta, blim.ml$eta)
dim(jac)     # 4095  583
qr(jac)$rank # 570  -> nicht identifizierbar ?!?

graded <- rbind(is.forward.graded(K), is.backward.graded(K))
rownames(graded) <- c("forward", "backward")

#           a  b  c  d  e  f  g  h  i  j   k   l
# item     10 18 31 38 62 65 68 70 84 89 117 127
# forward   1  0  1  1  0  1  1  0  1  0   0   0
# backward  0  1  1  0  1  1  0  0  0  1   1   1
