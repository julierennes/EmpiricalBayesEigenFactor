
# Importation of the citation matrix
rm(list = objects())
setwd("~/Data")
MatC <- read.csv("Matrix C Varin Cattelan Firth.csv", header  = T, sep =";", row.names = 1)
summary(MatC)
n <- nrow(MatC)
MatC <- t(MatC)

# Importation of the number of references published
References_ai <-  read.csv("ReferencesPublished", header  = T, sep =";")
summary(References_ai)
Teleportation_pi <-(References_ai/sum(References_ai)) # \pi = ai/a+

# Creation of the transition probability matrix P with pij=cij/ci+
MatP <-  sweep(MatC, 1, apply(MatC, 1, sum), FUN = "/")

###########
# EIFA
###########

MatC_withoutDiag <- MatC
diag(MatC_withoutDiag) <- 0
MatP_withoutDiag <-  sweep(MatC_withoutDiag, 1, apply(MatC_withoutDiag, 1, sum), FUN = "/")

alpha <- 0.85 # damping factor
MatG_EIFA <- alpha * MatP_withoutDiag + (1-alpha) * (rep(1, n)%*% t(Teleportation_pi))

MatG_EIFA_r  <- eigen(t(MatG_EIFA))$vector[,1]
MatG_EIFA_r  <- MatG_EIFA_r/sum(MatG_EIFA_r) 

# Equation 3 
MatG_EIFA_rstar <- (MatG_EIFA_r -(1-alpha)*(Teleportation_pi))/alpha 
resEIFA <- cbind.data.frame(rownames(MatC), as.numeric((MatG_EIFA_rstar/sum(MatG_EIFA_rstar) * 1000)[,1]))
colnames(resEIFA) <- c("Journal", "EIFA")
# Table 2 column 3, total influence score
res_EIFA <- resEIFA[order(resEIFA[,2], decreasing = T), ] 

res_EIFA_AI <-  cbind.data.frame(rownames(MatC), resEIFA[,2]/Teleportation_pi/1000)
colnames(res_EIFA_AI) <- c("Journal", "EIFA_AI")
# Table 3 column 3, article influence score
res_EIFA_AI <-res_EIFA_AI[order(res_EIFA[,2], decreasing = T), ]

#########
# EBEF
#########

# initialization \gamma_j^{0} = N c+j/c++ 
estim_gamma <- n*apply(MatC, 2, sum)/ sum(apply(MatC, 2, sum)) 
MatC_withoutdiag <- MatC
diag(MatC_withoutdiag) <- NA
# ni= \sum_{i \neq j} c_{ij}
ni <- apply(MatC_withoutdiag, 1, sum, na.rm = TRUE) 

## Fixed Point algorithm - Equation 15
NUMPART = DENOM = 0
 for ( l in 1:1000){  
  for ( j in 1:n){
  K = sum(estim_gamma)
  NUMPART = DENOM = 0
    for (i in 1:n){ 
     if (i!=j) {
       NUMPART = NUMPART + digamma(MatC_withoutdiag[i , j] + estim_gamma[j])
       DENOM =  DENOM  + digamma(ni[i] + (K - estim_gamma[i]))  - digamma(K - estim_gamma[i])   
                }
                  }
    NUM  = NUMPART  - (n-1) * digamma(estim_gamma[j])
    estim_gamma[j] <- estim_gamma[j] * NUM/DENOM
   }
}

RES_gamma <- cbind.data.frame(rownames(MatC),   estim_gamma)
colnames(RES_gamma) <- c("Journal", "gamma")
RES_gamma <- RES_gamma[as.numeric(rownames(res_EIFA)), ]

Kfinal <- sum(estim_gamma)
Kwithout_i <- Kfinal - estim_gamma
# Equation 8
Alpha_i <- ni/(ni+Kwithout_i)
A <- diag(Alpha_i)
mean(Alpha_i)

Teleportation_pistar <- sweep(matrix(estim_gamma,n,n, byrow = TRUE), 1, Kwithout_i, FUN = "/")
diag(Teleportation_pistar) <- 0

# Equation 9
MatG_EBEF <-t(MatP_withoutDiag) %*% A  + (t(Teleportation_pistar))%*% (diag(1, n)-A)

# First eigenvector of G^{\star}
MatG_EBEF_r <- eigen((MatG_EBEF))$vector[, 1]
resEBEF <- cbind.data.frame(rownames(MatC), as.numeric((MatG_EBEF_r /sum(MatG_EBEF_r ) * 1000)))
colnames(resEBEF) <- c("Journal", "EBEF")
# Table 2, column 4, total influence score
res_EBEF <- resEBEF[order(resEIFA[,2], decreasing = T), ]
res_EBEF

res_EBEF_AI <-  cbind.data.frame(rownames(MatC), resEBEF[,2]/Teleportation_pi/1000)
colnames(res_EBEF_AI) <- c("Journal", "EBEF_AI")
# Table 3 column 4, article influence score
res_EBEF_AI <-res_EBEF_AI[order(res_EIFA[,2], decreasing = T), ]

#########
# EBPR 
#########
estim_gamma <- rep(1, n)
ni <- apply(MatC, 1, sum, na.rm = TRUE) 
NUMPART = DENOM = 0
for ( l in 1:1000){  
for ( j in 1:n) {
  K = sum(estim_gamma)
  NUMPART = DENOM = 0
    for ( i in 1:n){ 
        NUMPART = NUMPART + digamma(MatC[i , j] + estim_gamma[j])
        DENOM =  DENOM  + digamma (ni[i] + (K - estim_gamma[i]))  - digamma(K - estim_gamma[i])   
    }
    NUM  = NUMPART  - (n-1) * digamma(estim_gamma[j])
    estim_gamma[j] <- estim_gamma[j] * NUM/DENOM
  }
}

Kfinal <- sum(estim_gamma)
Alpha_i <- ni/(ni+K)
A <- diag(Alpha_i)
Teleportation_pistar <- matrix(estim_gamma,n,n, byrow= TRUE)/ K


MatG_EBPR <- t(MatP)%*%A + (t(Teleportation_pistar))%*% (diag(1, n)-A)
MatG_EBPR_r <- eigen((MatG_EBPR))$vector[,1]

resEBPR <- cbind.data.frame(rownames(MatC), as.numeric((MatG_EBPR_r /sum(MatG_EBPR_r ) * 1000)))
colnames(resEBPR) <- c("Journal", "EBPR")
res_EBPR <- resEBPR[order(resEIFA[,2], decreasing = T), ]
res_EBPR

## Inversion method
estim_gamma <- n*  apply(MatC, 2, sum)/ sum(apply(MatC, 2, sum)) 
ni <- apply(MatC_withoutdiag, 1, sum, na.rm = TRUE) 
NUMPART = DENOM = 0
for ( l in 1:1000){
  for ( j in 1:n) {
    
    K = sum(estim_gamma)
    NUMPART = DENOM = 0
    for ( i in 1:n){ 
      if (i!=j) {  NUMPART = NUMPART + digamma(K - estim_gamma[i])  - digamma(ni[i] + (K - estim_gamma[i]))
      DENOM = DENOM +  digamma(MatC[i , j] + estim_gamma[j])
      }
    }  
    a = ((n-1)^(-1)) * (NUMPART + DENOM)
    estim_gamma[j] = estim_gamma[j] -  ( digamma(estim_gamma[j]) - a)/trigamma(estim_gamma[j])
  }
}
###

  