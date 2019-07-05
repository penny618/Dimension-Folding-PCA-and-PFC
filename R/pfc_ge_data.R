library(MASS)
library(compiler)
set.seed(2019)
source('C:/Users/Penny/Desktop/PFC/code/generate_data.R')

pL  <- pR <- 3
dL  <- dR <- 2
rR  <- rL <- 4
n   <- 50
tol <- 1e-10
omega <- matrix(c(0.50,-0.25,0,-0.25,0.5,-0.25,0,-0.25,0.5), 3, 3)
M <- matrix(c(0.886,0.266,0.062,0.266,0.248,0.048,0.062,0.048,0.015), 3, 3)

RawD <- pfcge.data(n, pL, pR, dL, dR, rL, rR, omega, M)  # dim(x)=pL*pR
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
omegahat <- diag(abs(rnorm(pR)))
Mhat <- diag(abs(rnorm(pL)))

# dimension folding-------------------------------------------------------------------------------------------
source('C:/Users/Penny/Desktop/code/pfc_ge.R')
source('C:/Users/Penny/Desktop/code/pfcge_exc.R')

n > max(pL/pR, pR/pL) - 1  # if true then Mhat and omegahat is invertible

if(n*pL > 30000 | n*pR > 30000) f <- 'pfcge.exc' else f <- 'pfc.ge'
f

# pfc.ge, dimension floding PFC with general error
pre <- pfc.ge(x, fy, gamma1, beta1, omegahat, tol)
tru <- list(RawD[4:7],omega=omega,M=M)  # true value
omegahat <- pre[['omegahat']]
Mhat     <- pre[['Mhat']]
gamma1   <- pre[['gamma1']]
gamma2   <- pre[['gamma2']]


# pfcge.exc, when npL*nPL or npR*npR are very large(>30000*30000)
pre <- pfcge.exc(x, fy, gamma1, beta1, omegahat, Mhat, tol)
tru <- list(RawD[4:7],omega=omega,M=M)  # true value
omegahat <- pre[['omegahat']]
Mhat     <- pre[['Mhat']]
gamma1   <- pre[['gamma1']]
gamma2   <- pre[['gamma2']]


# PCDF_Error-----------------------------------------------------------------------------
S1 <- kronecker(solve(omegahat)%*%gamma1, solve(Mhat)%*%gamma2)  # matrix(predict)
P1 <- S1 %*% solve(t(S1) %*% S1) %*% t(S1)  # projection matrix(predict)
S2 <- kronecker(solve(omega)%*%Gamma1, solve(M)%*%Gamma2)  # matrix
P2 <- S2 %*% solve(t(S2) %*% S2) %*% t(S2)  # projection matrix
P  <- P1 - P2
PCDF_Error <- (norm(P,type = 'F'))^2


# choice of dL and dR---------------------------------
BIC <- -2 * pre[['loglikelihood']] + log(n) * sum(c(pL*dL, dL*rL, dR*rR, pR*dR));BIC

AIC <- -2 * pre[['loglikelihood']] + 2 * sum(c(pL*dL, dL*rL, dR*rR, pR*dR));AIC








