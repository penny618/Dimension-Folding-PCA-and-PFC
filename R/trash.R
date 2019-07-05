setwd("C:/Users/Penny/Desktop/PFC/paper")

library(jpeg)
set.seed(123)
img <- readJPEG("C:/Users/Penny/Desktop/PFC/paper/data/dataset-resized/dataset-resized/glass/glass21.jpg")
#img.n <- readJPEG("C:/Users/Penny/Desktop/PFC/paper/data/dataset-resized/dataset-resized/glass/glass215.jpg",TRUE)
class(img)

imgRGB
plot(1:2,type = 'n')
rasterImage(img,1,1,2,2)
dim(img)
# 384 512   3
# 256  64 100
library(magick)
path <- "C:/Users/Penny/Desktop/PFC/paper/data/dataset-resized/dataset-resized/"
pathIn1 <- paste(path,"glass/glass",1:10,".jpg",sep = "")
pathIn2 <- paste(path,"metal/metal",1:10,".jpg",sep = "")
pathIn <- c(pathIn1, pathIn2)
n <- length(pathIn)

img <- image_read(pathIn)  # 读入图片
gray <- image_convert(img, colorspace='gray')  # 转为灰度图

gray_matrix <- array(0, c(384,512,n))
for(i in 1:10){
  gray_matrix[,,i] <- as.integer(image_data(gray[i]))
}
dim(gray_matrix)

n*pL*pR
pL <- 384
pR <- 512
dL <- dR <- 20
rR <- rL <- 1
#sig <- 1
tol <- 1e-4

x <- gray_matrix
dim(x)

y <- 10-1/2 
fy <- array(0, c(rL, rR, n))
for(i in 1:n){
  fy[,,i] <- matrix(y,1,1)
}


# dimension folding-------------------------------------------------------------------
covx <- matrix(rep(0, pR*pL), pR, pR)
for (i in n) {
  covx <- t(x[,,i] - apply(x, c(1, 2), mean)) %*% (x[,,i] - apply(x, c(1, 2), mean)) + covx
}
dim(covx)
# generate gamma1
gamma1 <- eigen(covx/n)$vectors[,1:dR]

# generate beta1
beta1  <- matrix(abs(rnorm(dR*rR)), dR, rR)

# generate omegahat
omegahat <- diag(abs(rnorm(pR)))

# generate Mhat
Mhat <- diag(abs(rnorm(pL)))

pre <- pfcge.exc(x, fy, gamma1, beta1, omegahat, Mhat, tol)
X11 <- pre[["X11"]]
dim(X11)
Y <- as.factor(c(rep("glass",10),rep("metal",10)))

# matrix straighten

X <- matrix(0, nrow = n, ncol = dL*dR)
for(i in 1:n){
  for(j in 1:dL){
    X[i,((j-1)*dR+1):(j*dR)] <- X11[j,,i]
  }
}
dim(X)

train <- data.frame(Y,X)

# plot
plot(density(X[11:20,1]),lty=1,col=1)#a
lines(density(X[1:10,1]),lty=2,col=2)#c
legend('topright',c('alcoholic group','control group'),lty = 1:2,col = 1:2)


# SVM
library(e1071)
svm.model <- svm(Y~.,data = train)
model.pred <- fitted(svm.model)
table(model.pred, Y)
model.pred

## pfcge.exc ############################################################################
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



save.image("picture.RData")







#####################################
library(magick)  # 加载magick包

image2chars <- function(pathIn='',
                        pathOutTxt=NULL,
                        pathOutImg=NULL,
                        jpg_quality=80,
                        width=100,
                        chars=c('&','#','w','s','k','d','t','j','i','.', ' '),
                        isColor=FALSE){
  ##### 参数
  # pathIn: 原始图片的路径，支持各种格式
  # pathOutTxt: 字符文本的输出路径，默认与原始图片在同一文件夹下，只是后缀为.txt;你也可指定其他路径
  # pathOutImg: 字符图片的输出路径，默认与原始图片在同一文件夹下，只是后缀为.jpg;你也可指定其他路径
  # jpg_quality: 字符图片的质量，范围0-100，越大图片越清晰，默认为80
  # width: 字符文本的宽度，默认为100，即一行100个字符；字符图片的尺寸也与其成正比
  # chars: 字符集，可自定义；默认为'&','#','w','s','k','d','t','j','i','.', ' '共11个字符
  # isColor: 字符图片是否为彩色，默认为黑白字符图片
  
  ##### 返回值
  # 无
  
  img <- image_read(pathIn)  # 读入图片
  gray <- image_convert(img, colorspace='gray')  # 转为灰度图
  rgb <- image_convert(img, colorspace='rgb')  # 转为rgb图
  
  ## 修改图片尺寸
  gray <- image_resize(gray, paste0(width,'x'))
  rgb <- image_resize(rgb, paste0(width,'x'))
  
  ## 获取图片灰度值矩阵，并将各像素值对应于相应字符
  gray <- as.integer(image_data(gray))[, , 1]
  w <- ncol(gray)   # 图片宽度
  h <- nrow(gray)  # 图片高度
  index <- findInterval(c(gray), seq(0, 255, length.out=length(chars)+1), rightmost.closed=T)
  labels <- chars[index]
  labels_mat <- matrix(labels, ncol=w)
  
  ## 输出字符文本，并保存成文件
  if(is.null(pathOutTxt))
    pathOutTxt <- paste0(pathIn,'.txt') # 文本文件名，与输入图片文件名一致，只是后缀为.txt
  write.table(labels_mat, pathOutTxt,
              quote=F, row.names=F,col.names=F)
  
  
  ## 绘制字符图片，给相应字符着色，并保存成文件
  if(isColor){
    rgb <- as.integer(image_data(rgb))
    r <- rgb[, , 1]  # red通道像素矩阵
    g <- rgb[, , 2]  # green通道像素矩阵
    b <- rgb[, , 3]  # blue通道像素矩阵
    
    cols <- rgb(c(r), c(g), c(b), maxColorValue=255) # 转化为颜色表示
  }
  
  if(is.null(pathOutImg))
    pathOutImg <- paste0(pathIn,'.jpg')  # 图片文件名，与输入图片文件名一致，只是后缀为.jpg
  jpeg(pathOutImg, width=16*w, height=16*h, quality=jpg_quality)
  op <- par(mar=c(0, 0, 0, 0))
  plot(0, 
       xlab='',
       ylab='',
       asp=1,
       xlim=c(0,w),
       ylim=c(0,h),
       xaxs="i",
       yaxs="i",
       type='n', 
       axes=FALSE)
  
  grid <- expand.grid(y=h:1-0.5, x=1:w-0.5)  # 各字符位置
  if(isColor){
    text(grid$x, grid$y, labels, cex=1.5, col=cols)  # 绘制彩色字符
  } else {
    text(grid$x, grid$y, labels, cex=1.5) # 绘制黑白字符
  }
  
  par(op)
  dev.off()
}

###################################






