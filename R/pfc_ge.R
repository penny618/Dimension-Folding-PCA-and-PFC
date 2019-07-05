# dimension floding PFC with general error

decomp <- function(mat, orders){
  decomp <- eigen(mat)$vectors %*% diag((eigen(mat)$values)^(orders)) %*% t(eigen(mat)$vectors)
  return(decomp)
}

pfc.ge <- function(x, fy, gamma1, beta1, omegahat, tol){
  XL <- matrix(0, n*pR, pL)
  FL <- matrix(0, n*pR, rL)
  XR <- matrix(0, n*pL, pR)
  FR <- matrix(0, n*pL, rR)
  L1 <- 100
  L0 <- 0
  
  while(abs(L1-L0) > tol){
    L0 <- L1
    # update gamma2 and beta2
    for(i in 1:n){
      XL[((i-1)*pR+1) : (i*pR),] <- t(x[,,i] %*% decomp(omegahat, -1/2))
      FL[((i-1)*pR+1) : (i*pR),] <- t(fy[,,i] %*% t(beta1) %*% t(gamma1) %*% decomp(omegahat, -1/2))
    } 
    #sigmaL <- t(XL) %*% FL%*%solve(t(FL)%*%FL)%*%t(FL) %*% XL/(n*pR)
    sigmaL <- t(XL) %*% FL %*% qr.solve(t(FL)%*%FL, t(FL), tol=1e-20) %*% XL/(n*pR)
    Mres   <- t(XL)%*%XL/(n*pR) - sigmaL
    UL     <- eigen(decomp(Mres,-1/2)%*%sigmaL%*%decomp(Mres,-1/2))$vectors
    lambdaL<- diag(eigen(decomp(Mres,-1/2)%*%sigmaL%*%decomp(Mres,-1/2))$values)
    DL     <- diag(c(rep(0,dL),(eigen(lambdaL)$values)[(dL+1):pL]))
    Mhat   <- Mres + decomp(Mres,1/2)%*%UL%*%DL%*%t(UL)%*%decomp(Mres,1/2) 
    
    gamma2 <- decomp(Mhat,1/2) %*% eigen(decomp(Mhat,-1/2)%*%sigmaL%*%decomp(Mhat,-1/2))$vectors[,1:dL]
    beta2  <- t(gamma2) %*% solve(Mhat) %*%t(XL) %*% FL %*% solve(t(FL)%*%FL)
    
    # update gamma1 and beta1
    for(i in 1:n){
      XR[((i-1)*pL+1) : (i*pL),] <- t(t(x[,,i]) %*% decomp(Mhat,-1/2))
      FR[((i-1)*pL+1) : (i*pL),] <- t(decomp(Mhat,-1/2)) %*% gamma2 %*% beta2 %*% fy[,,i]
    }
    #sigmaR  <- t(XR) %*% FR%*%solve(t(FR)%*%FR)%*%t(FR) %*% XR/(n*pL)
    sigmaR  <- t(XR) %*% FR %*% qr.solve(t(FR)%*%FR, t(FR), tol = 1e-20) %*% XR/(n*pL)
    omegares <- t(XR)%*%XR/(n*pL) - sigmaR
    UR <- eigen(decomp(omegares,-1/2)%*%sigmaR%*%decomp(omegares,-1/2))$vectors
    lambdaR <- diag(eigen(decomp(omegares,-1/2)%*%sigmaR%*%decomp(omegares,-1/2))$values)
    DR <- diag(c(rep(0,dR),(eigen(lambdaR)$values)[(dR+1):pR]))
    omegahat <- omegares + decomp(omegares,1/2)%*%UR%*%DR%*%t(UR)%*%decomp(omegares,1/2)
    
    gamma1 <- decomp(omegahat,1/2)%*%eigen(decomp(omegahat,-1/2)%*%sigmaR%*%decomp(omegahat,-1/2))$vectors[,1:dR]
    beta1  <- t(gamma1) %*% solve(omegahat) %*% t(XR) %*% FR %*% solve(t(FR)%*%FR)
    
    xbar <- apply(x, c(1,2), mean)
    # update log likelihood function
    l <- 0
    for(i in 1:n){
      delta <- x[,,i] - xbar - gamma2 %*% beta2 %*% fy[,,i] %*% t(beta1) %*% t(gamma1)
      l <- l + sum(diag(solve(omegahat) %*% t(delta) %*% solve(Mhat) %*% delta))
    }
    L1 <- n*pL*pR/2*log(2*pi) + n*pL/2*log(det(omegahat)) + n*pR/2*log(det(Mhat)) + 1/2*l
  }
  # return estimators
  X11 <- array(0, c(dL, dR, n))  # New predictor after dimension folding
  for(i in 1:n){
    X11[,,i] <- t(gamma2) %*% solve(Mhat) %*% x[,,i] %*% solve(omegahat) %*% gamma1
  }
  est <- list(gamma1=gamma1, gamma2=gamma2, beta1=beta1, beta2=beta2, 
              omegahat=omegahat, Mhat=Mhat, loglikelihood=-L1, X11 = X11)
  return(est)
}


