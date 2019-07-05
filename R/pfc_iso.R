# dimension floding PFC with isotropic error

pfc.iso <- function(x, fy, gamma1, beta1, tol){
  XL <- matrix(0, n*dR, pL)
  FL <- matrix(0, n*dR, rL)
  XR <- matrix(0, n*dL, pR)
  FR <- matrix(0, n*dL, rR)
  L1 <- 100
  L0 <- 0
  
  while(abs(L1-L0) > tol){
    L0 <- L1
    # update gamma2 and beta2
    for(i in 1:n){
      XL[((i-1)*dR+1) : (i*dR),] <- t(x[,,i] %*% gamma1)
      FL[((i-1)*dR+1) : (i*dR),] <- t(fy[,,i] %*% t(beta1))
    }
    sigmaL  <- t(XL) %*% FL%*%solve(t(FL)%*%FL)%*%t(FL) %*% XL/n
    gamma2 <- eigen(sigmaL)$vectors[,1:dL]
    beta2  <- t(gamma2) %*% t(XL) %*% FL %*% solve(t(FL)%*%FL)
    
    # update gamma1 and beta1
    for(i in 1:n){
      XR[((i-1)*dL+1) : (i*dL),] <- t(t(x[,,i]) %*% gamma2)
      FR[((i-1)*dL+1) : (i*dL),] <- beta2 %*% fy[,,i]
    }
    sigmaR  <- t(XR) %*% FR%*%solve(t(FR)%*%FR)%*%t(FR) %*% XR/n
    gamma1 <- eigen(sigmaR)$vectors[,1:dR]
    beta1  <- t(gamma1) %*% t(XR) %*% FR %*% solve(t(FR)%*%FR)
    
    #xbar <- apply(x, c(1,2), mean)
    # update log likelihood function
    l <- 0
    for(i in 1:n){
      delta <- x[,,i] - gamma2 %*% beta2 %*% fy[,,i] %*% t(beta1) %*% t(gamma1)
      l <- l + sum(diag(t(delta) %*% delta))
    }
    sigmaHat <- l/(n*pL*pR)
    L1 <- n*pL*pR/2*(log(2*pi*sigmaHat)) + 1/(2*sigmaHat)*l
  }
  # return estimators
  X11 <- array(0, c(dL, dR, n))  # New predictor after dimension folding
  for(i in 1:n){
    X11[,,i] <- t(gamma2) %*% x[,,i] %*% gamma1
  }
  est <- list(gamma1=gamma1, gamma2=gamma2, beta1=beta1, beta2=beta2, 
              loglikelihood=-L1, sigmaHat=l/(n*pL*pR), X11=X11)
  return(est)
}


