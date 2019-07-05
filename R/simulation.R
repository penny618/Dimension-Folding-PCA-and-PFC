library(MASS)
library(compiler)
source('C:/Users/Penny/Desktop/PFC/code/generate_data.R')
source('C:/Users/Penny/Desktop/PFC/code/pfc_ge.R')
source('C:/Users/Penny/Desktop/PFC/code/pfc_iso.R')

set.seed(2019)

# Figure2 PCDF Error(DF-PFC general error)-------------------------------------------------------------
pL  <- pR <- 3
dL  <- dR <- 2
rR  <- rL <- 4
N   <- c(100, 150)
#N   <- c(30, 50, 80, 100, 150)
tol <- 1e-10
omega <- matrix(c(0.50,-0.25,0,-0.25,0.5,-0.25,0,-0.25,0.5), 3, 3)
M <- matrix(c(0.886,0.266,0.062,0.266,0.248,0.048,0.062,0.048,0.015), 3, 3)

Times <- 200
Nc <- length(N)
Nr <- Times/length(N)
PCDF_Error.ge2 <- matrix(rep(0, Nr*Nc), nrow = Nr, ncol = Nc)
colnames(PCDF_Error.ge2) <- N

for(j in 1:Times){
  nc <- (j-1) %/% Nr + 1  # when n = N[ith]
  nr <- j - (nc - 1) * 100   # the j th repeat
  n  <- N[nc]
  
  RawD <- pfcge.data(n, pL, pR, dL, dR, rL, rR, omega, M)  # dim(x)=pL*pR
  x      <- RawD[['x']]
  fy     <- RawD[['fy']]
  Gamma1 <- RawD[['Gamma1']]
  Gamma2 <- RawD[['Gamma2']]
  # Generate gamma1 and beta1 
  covx <- matrix(rep(0, pR*pL), pR, pR)
  for (i in n) {
    covx <-  t(x[,,i]) %*% x[,,i] + covx
  }
  
  gamma1 <- eigen(covx/n)$vectors[,1:dR]
  beta1  <- matrix(abs(rnorm(dR*rR)), dR, rR)
  omegahat <- diag(abs(rnorm(pR)))
  Mhat <- diag(abs(rnorm(pL)))
  
  # pfc.ge, dimension floding PFC with general error
  pre <- pfc.ge(x, fy, gamma1, beta1, omegahat, tol)
  tru <- list(RawD[4:7], omega = omega, M = M)  # true value
  omegahat <- pre[['omegahat']]
  Mhat     <- pre[['Mhat']]
  gamma1   <- pre[['gamma1']]
  gamma2   <- pre[['gamma2']]
  
  
  # PCDF_Error
  S1 <- kronecker(solve(omegahat)%*%gamma1, solve(Mhat)%*%gamma2)  # matrix(predict)
  P1 <- S1 %*% solve(t(S1) %*% S1) %*% t(S1)  # projection matrix(predict)
  S2 <- kronecker(solve(omega)%*%Gamma1, solve(M)%*%Gamma2)  # matrix
  P2 <- S2 %*% solve(t(S2) %*% S2) %*% t(S2)  # projection matrix
  P  <-  P1 - P2
  PCDF_Error.ge2[nr,nc] <- (norm(P, type = 'F'))^2
  
  
  print(paste("Come on, I just finished the", j, "th replication."))

}
setwd('C:/Users/Penny/Desktop/code')
save.image("error_ge2.RData")
pdf(file = '../code/fig/error_ge2.pdf')
boxplot(PCDF_Error.ge2)
dev.off()




# Figure3 Estimation error(isotropic error)-------------------------------------------------------------
pL  <- pR <- 10
dL  <- dR <- 2
rR  <- rL <- 4
N   <- c(120, 150, 200, 300, 500)
tol <- 1e-10
sig <- 0.8

Times <- 6
Nc <- length(N)
Nr <- Times/length(N)
PCDF_Error.iso <- matrix(rep(0, Nr*Nc), nrow = Nr, ncol = Nc)
colnames(PCDF_Error.iso) <- N
for(j in 1:Times){
  nc <- (j-1) %/% Nr + 1  # when n = N[ith]
  nr <- j - (nc - 1) * 100  # the j th repeat
  n  <- N[nc]
  
  RawD <- pfciso.data(n,pL,pR,dL,dR,rL,rR)  # dim(x)=pL*pR
  x      <- RawD[['x']]
  fy     <- RawD[['fy']]
  Gamma1 <- RawD[['Gamma1']]
  Gamma2 <- RawD[['Gamma2']]
  # Generate gamma1 and beta1 
  covx <- matrix(rep(0, pR*pL), pR, pR)
  for (i in n) {
    covx <- t(x[,,i]) %*% x[,,i] + covx
  }
  
  gamma1 <- eigen(covx/n)$vectors[,1:dR]
  beta1  <- matrix(abs(rnorm(dR*rR)), dR, rR)
  
  # pfc.iso, dimension floding PFC with isotropic error
  pre <- pfc.iso(x, fy, gamma1, beta1, tol)
  tru <- RawD[4:7]
  gamma1   <- pre[['gamma1']]
  gamma2   <- pre[['gamma2']]
  
  # PCDF_Error
  S1 <- kronecker(gamma1, gamma2)  # matrix(predict)
  P1 <- S1 %*% solve(t(S1) %*% S1) %*% t(S1)  # projection matrix(predict)
  S2 <- kronecker(Gamma1, Gamma2)  # matrix
  P2 <- S2 %*% solve(t(S2) %*% S2) %*% t(S2)  # projection matrix
  P  <- P1 - P2
  
  PCDF_Error.iso[nr,nc] <- (norm(P,type = 'F'))^2
  
  
  print(paste("Come on, I just finished the", j, "th replication."))
  
}

setwd('C:/Users/Penny/Desktop/code')
save.image("error_iso.RData")
pdf(file = '../code/fig/error_iso.pdf')
boxplot(PCDF_Error.iso)
dev.off()


# check(isotropic error)-------------------------------------------------------------
source('C:/Users/Penny/Desktop/PFC/code/generate_data.R')
source('C:/Users/Penny/Desktop/PFC/code/pfc_ge.R')
source('C:/Users/Penny/Desktop/PFC/code/pfc_iso.R')

p <- c(5,10,15,20)
dL  <- dR <- 2
rR  <- rL <- 4
n   <- 100
sig <- 0.8
tol <- 1e-10

Times <- 400
Nc <- length(p)
Nr <- Times/length(p)
PCDF_Error.iso <- matrix(rep(0, Nr*Nc), nrow = Nr, ncol = Nc)
colnames(PCDF_Error.iso) <- p
for(j in 1:Times){
  nc <- (j-1) %/% Nr + 1  # when n = N[ith]
  nr <- j - (nc - 1) * 100  # the j th repeat
  pL  <- pR <- p[nc]
  
  
  RawD <- pfciso.data(n,pL,pR,dL,dR,rL,rR)  # dim(x)=pL*pR
  x      <- RawD[['x']]
  fy     <- RawD[['fy']]
  Gamma1 <- RawD[['Gamma1']]
  Gamma2 <- RawD[['Gamma2']]
  # Generate gamma1 and beta1 
  covx <- matrix(rep(0, pR*pL), pR, pR)
  for (i in n) {
    covx <- t(x[,,i]) %*% x[,,i] + covx
  }
  
  gamma1 <- eigen(covx/n)$vectors[,1:dR]
  beta1  <- matrix(abs(rnorm(dR*rR)), dR, rR)
  
  # pfc.iso, dimension floding PFC with isotropic error
  pre <- pfc.iso(x, fy, gamma1, beta1, tol)
  gamma1   <- pre[['gamma1']]
  gamma2   <- pre[['gamma2']]
  
  # PCDF_Error
  S1 <- kronecker(gamma1, gamma2)  # matrix(predict)
  P1 <- S1 %*% solve(t(S1) %*% S1) %*% t(S1)  # projection matrix(predict)
  S2 <- kronecker(Gamma1, Gamma2)  # matrix
  P2 <- S2 %*% solve(t(S2) %*% S2) %*% t(S2)  # projection matrix
  P  <- P1 - P2
  
  PCDF_Error.iso[nr,nc] <- (norm(P,type = 'F'))^2
  
  
  print(paste("Come on, I just finished the", j, "th replication."))
  
}

setwd('C:/Users/Penny/Desktop/PFC/code')
save.image("error_iso.RData")
pdf(file = '../code/fig/error_iso_check.pdf')
boxplot(PCDF_Error.iso)
dev.off()
PCDF_Error.iso
colMeans(PCDF_Error.iso)
summary(PCDF_Error.iso)
