library(compiler)
# pfciso----------------------------------------------------------------------------
pfciso.data <- function(n, pL, pR, dL, dR, rL, rR){
  Gamma1 <- matrix(rnorm(pR*dR, 0, 1), pR, dR)
  Gamma2 <- matrix(rnorm(pL*dL, 0, 1), pL, dL)
  
  y <- rnorm(n)  
  fy <- array(0, c(rL, rR, n))
  for(i in 1:n){
    fy[,,i] <- diag(c(y[i], (y[i])^2, (y[i])^3, (y[i])^4))
  }
  
  Beta1 <- matrix(rnorm(dR*rR, 1, 2), dR, rR)
  Beta2 <- matrix(abs(rnorm(dL*rL, 2, 2)), dL, rL)
  
  x <- array(0, c(pL, pR, n))
  
  e <- mvrnorm(n, rep(0, pL*pR), diag(rep(1, pL)) %x% diag(rep(1, pR)))
  e <- t(e)
  for(i in 1:n){
    x[,,i] <- Gamma2 %*% Beta2 %*% fy[,,i] %*% t(Beta1) %*% t(Gamma1) + matrix(e[,i], pL, pR)
  }
  for(i in 1:n){                   
    x[,,i] <- x[,,i] - apply(x, c(1, 2), mean)  
  }
  RawD <- list(x=x,y=y,fy=fy,Gamma1=Gamma1,Gamma2=Gamma2,Beta1=Beta1,Beta2=Beta2)
  return(RawD)
}
pfciso.data <- cmpfun(pfciso.data)
# data generation with general error 
pfcge.data <- function(n, pL, pR, dL, dR, rL, rR, omega, M){
  Gamma1 <- matrix(rnorm(pR*dR, 0, 1), pR, dR)
  Gamma2 <- matrix(rnorm(pL*dL, 0, 1), pL, dL)
  
  y <- rnorm(n)  
  fy <- array(0, c(rL, rR, n))
  for(i in 1:n){
    fy[,,i] <- diag(c(y[i], (y[i])^2, (y[i])^3, (y[i])^4))
  }
  
  Beta1 <- matrix(rnorm(dR*rR, 1, 2), dR, rR)
  Beta2 <- matrix(abs(rnorm(dL*rL, 2, 2)), dL, rL)
  
  x <- array(0, c(pL, pR, n))
  e <- mvrnorm(n, rep(0, pL*pR), omega %x% M)
  e <- t(e)
  for(i in 1:n){
    x[,,i] <- Gamma2 %*% Beta2 %*% fy[,,i] %*% t(Beta1) %*% t(Gamma1) + matrix(e[,i], pL, pR)
  }
  for(i in 1:n){                   
    x[,,i] <- x[,,i] - apply(x, c(1, 2), mean)  
  }
  RawD <- list(x=x,y=y,fy=fy,Gamma1=Gamma1,Gamma2=Gamma2,Beta1=Beta1,Beta2=Beta2)
  return(RawD)
  
}
pfcge.data <- cmpfun(pfcge.data)









