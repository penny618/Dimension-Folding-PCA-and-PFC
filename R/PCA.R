library(MASS)
set.seed(2019)

# data generation ----------------------------------------------------------------------------
pca.data <- function(n, pL, pR, dL, dR, sig){
  Gamma1 <- matrix(rnorm(pR*dR, 0, 1), pR, dR)
  Gamma2 <- matrix(rnorm(pL*dL, 0, 1), pL, dL)
  y <- array(rnorm(dL*dR, 1, 2) , c(dL,dR,n))
  x <- array(0, c(pL, pR, n))
  e <- mvrnorm(n, rep(0, pL*pR), sig * (diag(rep(1, pL)) %x% diag(rep(1, pR))))
  e <- t(e)
  for(i in 1:n){
    x[,,i] <- Gamma2 %*% y[,,i] %*%  t(Gamma1) + matrix(e[,i], pL, pR)
  }
  for(i in 1:n){                   
    x[,,i] <- x[,,i] - apply(x, c(1, 2), mean)  
  }
  for(i in 1:n){                   
    y[,,i] <- y[,,i] - apply(y, c(1, 2), mean)  
  }
  RawD <- list(x=x,y=y, Gamma1=Gamma1,Gamma2=Gamma2 )
  return(RawD)
}

p <- c(5, 10, 15, 20, 30)
pL  <- pR <- p[1]
dL  <- dR <- 2
sig <- 1
n   <- 100
tol <- 1e-10


RawD <- pca.data(n,pL,pR,dL,dR,sig)  # dim(x)=pL*pR
x <- RawD[['x']]
y <- RawD[['y']]
Gamma1 <- RawD[['Gamma1']]
Gamma2 <- RawD[['Gamma2']]

# Generate gamma1  
covx <- matrix(rep(0, pR*pL), pR, pR)
for (i in n) {
  covx <- t(x[,,i]) %*% x[,,i] + covx
}
# dimension floding PCA with isotropic error 
pca.iso <- function(x, y, gamma1, tol){
  L1 <- 100
  L0 <- 0
  SigmaL <- matrix(0,pL,pL)
  SigmaR <- matrix(0,pR,pR)
  
  while(abs(L1-L0) > tol){
    L0 <- L1
    # update gamma2 
    for (i in 1:n){
      sigmaL <- x[,,i]%*% gamma1%*%t(gamma1)%*%t(x[,,i])
      SigmaL <- SigmaL+sigmaL
    }
    SigmaL <- SigmaL/n
    gamma2 <- eigen(SigmaL)$vectors[,1:dL]
    # update gamma1 
    for (i in 1:n){
      sigmaR<- t(x[,,i])%*% gamma2%*%t(gamma2)%*%x[,,i]
      SigmaR <- SigmaR+sigmaR
    }
    SigmaR <- SigmaR/n
    gamma1 <- eigen(SigmaR)$vectors[,1:dR]
    # update log likelihood function
    delta <- x[,,1] - gamma2 %*%y[,,1] %*%t(gamma1)
    l <- sum(diag(t(delta) %*% delta))
    for(i in 2:n){
      delta <- x[,,i] - gamma2 %*%y[,,i] %*%t(gamma1)
      l <- l + sum(diag(t(delta) %*% delta))
    }
    sigmaHat <- l/(n*pL*pR)
    L1 <- n*pL*pR/2*(log(2*pi*sigmaHat)) + 1/(2*sigmaHat)*l
  }
  # return estimators
  est <- list(gamma1=gamma1, gamma2=gamma2, loglikelihood=-L1,sigmaHat=sigmaHat)
  return(est)
}

# dimension folding PCA-------------------------------------------------------------------------
gamma1 <- eigen(covx/n)$vectors[,1:dR]
pre <- pca.iso(x, y, gamma1,  tol)
tru <- RawD[3:4]
gamma1   <- pre[['gamma1']]
gamma2   <- pre[['gamma2']]

# PCDF_Error-----------------------------------------------------------------------------
S1 <- kronecker(gamma1, gamma2)  # matrix(predict)
P1 <- S1 %*% solve(t(S1) %*% S1) %*% t(S1)  # projection matrix(predict)
S2 <- kronecker(Gamma1, Gamma2)  # matrix
P2 <- S2 %*% solve(t(S2) %*% S2) %*% t(S2)  # projection matrix
P  <- P1 - P2
PCDF_Error <- (norm(P,"F"))^2
PCDF_Error

# simulation----------------------------------------------------------------------------
p <- c(5, 10, 15, 20, 30)

dL  <- dR <- 2
sig <- 1
n   <- 100
tol <- 1e-10

Times <- 500
Nc <- length(p)
Nr <- Times/length(p)
PCDF_Error.pca <- matrix(rep(0, Nr*Nc), nrow = Nr, ncol = Nc)
colnames(PCDF_Error.pca) <- p

for(j in 1:Times){
  nc <- (j-1) %/% Nr + 1  # when n = N[ith]
  nr <- j - (nc - 1) * 100  # the j th repeat
  pL  <- pR <- p[nc]
  
  RawD <- pca.data(n,pL,pR,dL,dR,sig)  # dim(x)=pL*pR
  x <- RawD[['x']]
  y <- RawD[['y']]
  Gamma1 <- RawD[['Gamma1']]
  Gamma2 <- RawD[['Gamma2']]
  # Generate gamma1  
  covx <- matrix(rep(0, pR*pL), pR, pR)
  for (i in n) {
    covx <- t(x[,,i]) %*% x[,,i] + covx
  }
  
  # pca.iso, dimension floding PCA with isotropic error
  gamma1 <- eigen(covx/n)$vectors[,1:dR]
  pre <- pca.iso(x, y, gamma1,  tol)
  tru <- RawD[3:4]
  gamma1   <- pre[['gamma1']]
  gamma2   <- pre[['gamma2']]
  
  # PCDF_Error
  S1 <- kronecker(gamma1, gamma2)  # matrix(predict)
  P1 <- S1 %*% solve(t(S1) %*% S1) %*% t(S1)  # projection matrix(predict)
  S2 <- kronecker(Gamma1, Gamma2)  # matrix
  P2 <- S2 %*% solve(t(S2) %*% S2) %*% t(S2)  # projection matrix
  P  <- P1 - P2
  PCDF_Error.pca[nr,nc] <- (norm(P, type = 'F'))^2
  
  
  print(paste("Come on, I just finished the", j, "th replication."))
  
}
setwd('C:/Users/Penny/Desktop/code')
save.image("error_pca.RData")
pdf(file = '../code/fig/error_pca.pdf')
boxplot(PCDF_Error.pca)
dev.off()





