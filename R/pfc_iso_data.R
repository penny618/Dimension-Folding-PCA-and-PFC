library(MASS)
set.seed(2019)

source('C:/Users/Penny/Desktop/PFC/code/R/generate_data.R')

pL  <- pR <- 10
dL  <- dR <- 2
rR  <- rL <- 4
sig <- 0.8
n   <- 120
tol <- 1e-10

RawD <- pfciso.data(n,pL,pR,dL,dR,rL,rR)  # dim(x)=pL*pR
x      <- RawD[['x']]
fy     <- RawD[['fy']]
Gamma1 <- RawD[['Gamma1']]
Gamma2 <- RawD[['Gamma2']]
# Generate gamma1 and beta1 ------------------------------------------------------------------
covx <- matrix(rep(0, pR*pL), pR, pR)
for (i in n) {
  covx <- t(x[,,i]) %*% x[,,i] + covx
}

gamma1 <- eigen(covx/n)$vectors[,1:dR]
beta1  <- matrix(abs(rnorm(dR*rR)), dR, rR)
# dimension folding-------------------------------------------------------------------------------------------
source('C:/Users/Penny/Desktop/PFC/code/pfc_iso.R')
pre <- pfc.iso(x, fy, gamma1, beta1, tol)
tru <- RawD[4:7]
gamma1   <- pre[['gamma1']]
gamma2   <- pre[['gamma2']]

# PCDF_Error-----------------------------------------------------------------------------
S1 <- kronecker(gamma1, gamma2)  # matrix(predict)
P1 <- S1 %*% solve(t(S1) %*% S1) %*% t(S1)  # projection matrix(predict)
S2 <- kronecker(Gamma1, Gamma2)  # matrix
P2 <- S2 %*% solve(t(S2) %*% S2) %*% t(S2)  # projection matrix
P  <- P1 - P2

PCDF_Error <- (norm(P,type = 'F'))^2
