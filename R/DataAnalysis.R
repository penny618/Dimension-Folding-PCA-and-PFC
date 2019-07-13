library(eegkitdata)
data("eegdata")
str(eegdata)

pL <- 256
pR <- 64
dL <- dR <- 1
rR <- rL <- 1
n  <- 100
tol <- 1e-5

table(eegdata$time)
time <- eegdata$time + 1
#table(time)

channel <- eegdata$channel
Y <- eegdata$group
table(Y)/sum(table(Y))
y <- matrix(as.numeric(table(Y)[1]-(table(Y)/sum(table(Y)))[1]),1,1); y # dim(fy)=1*1
fy <- array(0, c(rL, rR, n))
for(i in 1:n){
  fy[,,i] <- y
}

names_channel <- eegdata$channel[256 * ((1:n) - 1) + 1]
names_channel

x <- array(0, c(pL, pR, n))
for(l in 1:n){
  eeg <- eegdata[(256*64*(l-1)+1):(256*64*l),]
  for(j in 1:64){
    x[,j,l] <- eeg[(256*(j-1)+1):(256*j),7]
  }
}
dim(x)  # 256  64 100

source('C:/Users/Penny/Desktop/PFC/code/R/pfcge_exc.R')
source('C:/Users/Penny/Desktop/PFC/code/R/pfc_ge.R')

# Generate gamma1 and beta1 ------------------------------------------------------------------
covx <- matrix(rep(0, pR*pL), pR, pR)
for (i in n) {
  covx <- t(x[,,i] - apply(x, c(1, 2), mean)) %*% (x[,,i] - apply(x, c(1, 2), mean)) + covx
}
dim(covx)

gamma1 <- eigen(covx/n)$vectors[,1:dR]
beta1  <- matrix(abs(rnorm(dR*rR)), dR, rR)
omegahat <- diag(abs(rnorm(pR)))
Mhat <- diag(abs(rnorm(pL)))
# pfcge.exc, when npL*nPL or npR*npR are very large(>30000*30000)
pre <- pfc.ge(x, fy, gamma1, beta1, omegahat, tol)




pre <- pfcge.exc(x, fy, gamma1, beta1, omegahat, Mhat, tol)
omegahat <- pre[['omegahat']]
Mhat     <- pre[['Mhat']]
gamma1   <- pre[['gamma1']]
gamma2   <- pre[['gamma2']]
X11      <- pre[['X11']]

xx <- rep(0,100)
for(i in 1:100){
  xx[i] <- X11[1,1,i]
}
density(xx)
plot(density(xx[1:50]),lty=1,col=1)#a
lines(density(xx[51:100]),lty=2,col=2)#c
legend('topright',c('alcoholic group','control group'),lty = 1:2,col = 1:2)

save.image('eeg.RData')
# quadratic discriminant analysis-----------------------------------------
library(MASS)
library(ggplot2)
model1=lda(Species~.,data=iris,prior=c(1,1,1)/3)
table(Species,predict(model1)$class)

ld=predict(model1)$x
p=ggplot(cbind(iris,as.data.frame(ld))
         ,aes(x=LD1,y=LD2))
p+geom_point(aes(colour=Species),alpha=0.8,size=3)


model2=qda(Species~.,data=iris,cv=T)
predict(model)$posterior

yy <- rep(NA,100)
yy[1:50] <- 'a'
yy[51:100] <- 'c'
table(yy)

qeeg <- data.frame(yy,xx)
str(qeeg)
qda.fit <- qda(yy ~ ., data = qeeg, cv = T)
qda.pred<-predict(qda.fit,qeeg)
qda.pred$class
qda.pred$posterior[1:10,]
table(qda.pred$class,qeeg$yy)

122-107  # paper correctly classified 107 subjects out of the total 122 subjects
(122-107)/122
(122-94)/122

(26+7)/100

## pfc.ge #############################################################################
# dimension floding PFC with general error

decomp <- function(mat, orders){
  mat <- (mat+t(mat))/2
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
    L1 <- l
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

## pfcge.exc ############################################################################
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
    L1 <- l
    
  }
  # return estimators
  X11 <- array(0, c(dL, dR, n))
  for(i in 1:n){
    X11[,,i] <- t(gamma2) %*% solve(Mhat) %*% x[,,i] %*% solve(omegahat) %*% gamma1
  }
  est <- list(gamma1=gamma1, gamma2=gamma2, beta1=beta1, beta2=beta2, 
              omegahat=omegahat, Mhat=Mhat,loglikelihood=-L1, X11=X11)
  return(est)
}



