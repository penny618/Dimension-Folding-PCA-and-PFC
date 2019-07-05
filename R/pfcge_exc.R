# when npL*nPL or npR*npR are very large(>30000*30000) and exceed the storage limit in R software.

decomp <- function(mat, orders){
  decomp <- eigen(mat)$vectors %*% diag((eigen(mat)$values)^(orders)) %*% t(eigen(mat)$vectors)
  return(decomp)
}

pfcge.exc <- function(x, fy, gamma1, beta1, omegahat, Mhat, tol){
  XL <- matrix(0, n*pR, pL)
  FL <- matrix(0, n*pR, rL)
  XR <- matrix(0, n*pL, pR)
  FR <- matrix(0, n*pL, rR)
  z  <- array(0, c(pL, pR, n))
  
  L1 <- 100
  L0 <- 0
  
  while(abs(L1-L0) > tol){
    L0 <- L1
    # standardizes the predictors
    for(i in 1:n){
      z[,,i] <- decomp(Mhat,-1/2) %*% x[,,i] %*% decomp(omegahat,-1/2)
    }
    
    # update gamma2 and beta2
    gamma1z <- decomp(omegahat,-1/2) %*% gamma1
    for(i in 1:n){
      XL[((i-1)*dR+1) : (i*dR),] <- t(z[,,i] %*% gamma1z)
      FL[((i-1)*dR+1) : (i*dR),] <- t(fy[,,i] %*% t(beta1))
    }
    sigmaL  <- t(XL) %*% FL%*%solve(t(FL)%*%FL)%*%t(FL) %*% XL/n
    gamma2z <- eigen(sigmaL)$vectors[,1:dL]
    gamma2 <- decomp(Mhat,1/2) %*% gamma2z
    beta2  <- t(gamma2z) %*% t(XL) %*% FL %*% solve(t(FL)%*%FL)
    
    # update gamma1 and beta1
    for(i in 1:n){
      XR[((i-1)*dL+1) : (i*dL),] <- t(t(z[,,i]) %*% gamma2z)
      FR[((i-1)*dL+1) : (i*dL),] <- beta2 %*% fy[,,i]
    }
    sigmaR  <- t(XR) %*% FR%*%solve(t(FR)%*%FR)%*%t(FR) %*% XR/n
    gamma1z <- eigen(sigmaR)$vectors[,1:dR]
    gamma1 <- decomp(omegahat,1/2) %*% gamma1z
    beta1  <- t(gamma1z) %*% t(XR) %*% FR %*% solve(t(FR)%*%FR)
    
    xbar <- apply(x, c(1,2), mean)
    
    # update omegahat
    for(i in 1:n){
      delta <- t(x[,,i]-xbar-gamma2 %*% beta2 %*% fy[,,i] %*% t(beta1) %*% t(gamma1)) %*% solve(Mhat) %*% (x[,,i]-xbar-gamma2 %*% beta2 %*% fy[,,i] %*% t(beta1) %*% t(gamma1))
      omegahat <- omegahat + delta
    }
    omegahat <- omegahat/(n*pL)
    
    # update Mhat
    for(i in 1:n){
      delta <- (x[,,i]-xbar-gamma2 %*% beta2 %*% fy[,,i] %*% t(beta1) %*% t(gamma1)) %*% solve(omegahat) %*% t(x[,,i]-xbar-gamma2 %*% beta2 %*% fy[,,i] %*% t(beta1) %*% t(gamma1))
      Mhat <- Mhat + delta
    }
    Mhat <- Mhat/(n*pR)
    
    # update log likelihood function
    l <- 0
    for(i in 1:n){
      delta <- x[,,i] - xbar - gamma2 %*% beta2 %*% fy[,,i] %*% t(beta1) %*% t(gamma1)
      l <- l + sum(diag(solve(omegahat) %*% t(delta) %*% solve(Mhat) %*% delta))
    }
    #L1 <- n*pL*pR/2*log(2*pi) + n*pL/2*log(det(omegahat)) + n*pR/2*log(det(Mhat)) + 1/2*l
    L1 <- l
    
  }
  # return estimators
  X11 <- array(0, c(dL, dR, n))
  for(i in 1:n){
    X11[,,i] <- t(gamma2) %*% solve(Mhat) %*% x[,,i] %*% solve(omegahat) %*% gamma1
  }
  est <- list(gamma1=gamma1, gamma2=gamma2, beta1=beta1, beta2=beta2, 
              omegahat=omegahat, Mhat=Mhat,loglikelihood=-(n*pL*pR/2*log(2*pi) + n*pL/2*log(det(omegahat)) + n*pR/2*log(det(Mhat)) + 1/2*L1), X11=X11)
  return(est)
}



