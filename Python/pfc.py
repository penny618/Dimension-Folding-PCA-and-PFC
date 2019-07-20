# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 12:42:29 2019

@author: Penny
"""



# In[] 获取当前工作路径
import os
os.getcwd()

# In[] 切换工作目录
os.chdir('C:/Users/Penny/Desktop/PFC/paper/data/dataset-resized/dataset-resized/')
os.getcwd()

# In[]
import numpy as np
import pandas as pd
import random as r

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

r.seed(2019)
# In[] 奇异值分解 A=UDV'
def decomp(mat, orders):
    mat = mat + mat.T
    mat_val, mat_vec = np.linalg.eig(mat)
    dec = np.dot(np.dot(mat_vec, np.identity(len(mat_val))*(mat_val**orders)), mat_vec.T)
    return dec

# In[]

# In[]  read images
nk=10
array_of_img = []
for i in ["glass/glass","paper/paper","metal/metal"]:
    for j in (np.arange(nk)+1):
        img = mpimg.imread(i + str(j) + ".jpg")
        img_grey = np.dot(img[...,:3], [0.299, 0.587, 0.114])  # convert rgb to grey
        array_of_img.append(img_grey)

n = len(array_of_img)
#plt.imshow(img)

# In[]
pL = array_of_img[0].shape[0]
pR = array_of_img[0].shape[1]
dL = dR = 20
rR = rL = 2
#sig = 1
tol = 1e-4
n*pL*pR
# In[] compute covx
covx = 0
for i in np.arange(n):
    covx += array_of_img[i].T.dot(array_of_img[i])

gamma1 = np.linalg.eig(covx)[1][:,:dR]  # generate gamma1
beta1 = np.array([r.random() for i in np.arange(dR*rR)]).reshape((dR,rR))  # generate beta1
omegahat = np.identity(pR)*[np.abs(r.random()) for i in np.arange(pR)]  # generate omegahat
Mhat = np.identity(pL)*[np.abs(r.random()) for i in np.arange(pL)]  # generate Mhat
# In[]
x = array_of_img
# standardizes the predictors
z = [decomp(Mhat,-1/2).dot(x[i]).dot(decomp(omegahat,-1/2)) for i in np.arange(n)]

y = nk - nk/n
fy = [np.eye(2)*[y,y**2] for i in np.arange(n)]

# In[]
XL = np.zeros((n*dR,pL))
FL = np.zeros((n*dR,rL))
XR = np.zeros((n*dL,pR))
FR = np.zeros((n*dL,rR))
L1 = 100
L0 = 0

while np.abs(L1-L0) > tol:
    delta = 9
    L0 = L1
    
    # standardizes the predictors
    z = [decomp(Mhat,-1/2).dot(x[i]).dot(decomp(omegahat,-1/2)) for i in np.arange(n)]
    
    # update gamma2 and beta2
    gamma1z = decomp(omegahat,-1/2).dot(gamma1)
    for i in (np.arange(n)+1):
        XL[((i-1)*dR) : (i*dR),:] = (z[i-1].dot(gamma1z)).T
        FL[((i-1)*dR) : (i*dR),:] = (fy[i-1].dot(beta1.T)).T

    sigmaL = (XL.T.dot(FL)).dot(np.mat(FL.T.dot(FL)).I).dot(FL.T).dot(XL)/n
    gamma2z = (np.linalg.eig(sigmaL+np.eye(pL)*delta)[1][:,:dL]).real
    #gamma2z = (np.linalg.eig(sigmaL)[1][:,:dL]).real
    gamma2 = decomp(Mhat,1/2).dot(gamma2z)
    beta2 = gamma2z.T.dot(XL.T).dot(FL).dot(np.mat(FL.T.dot(FL)).I)
    
    # update gamma1 and beta1
    for i in (np.arange(n)+1):
        XR[((i-1)*dL) : (i*dL),:] = (z[i-1].T.dot(gamma2z)).T
        FR[((i-1)*dL) : (i*dL),:] = beta2.dot(fy[i-1])
    sigmaR = (XR.T.dot(FR)).dot(np.mat(FR.T.dot(FR)).I).dot(FR.T).dot(XR)/n
    gamma1z = (np.linalg.eig(sigmaR+np.eye(pR)*delta)[1][:,:dR]).real
    #gamma1z = (np.linalg.eig(sigmaR)[1][:,:dR]).real
    gamma1 = decomp(omegahat,1/2).dot(gamma1z)
    beta1 = gamma1z.T.dot(XR.T).dot(FR).dot(np.mat(FR.T.dot(FR)).I)

    xbar = np.apply_over_axes(np.mean,np.array(x),0)
    
    # update omegahat
    for i in (np.arange(n)+1):
        delta = (x[i-1] - xbar).reshape(pL,pR) - gamma2.dot(beta2).dot(fy[i-1]).dot(beta1.T).dot(gamma1.T)
        temp = delta.T.dot(np.mat(Mhat).I).dot(delta)
        omegahat += temp 
    omegahat = omegahat/(n*pL)
    
    # update Mhat
    for i in (np.arange(n)+1):
        delta = (x[i-1] - xbar).reshape(pL,pR) - gamma2.dot(beta2).dot(fy[i-1]).dot(beta1.T).dot(gamma1.T)
        temp = delta.dot(np.mat(omegahat).I).dot(delta.T)
        Mhat += temp 
    Mhat = Mhat/(n*pR)
    
    # update log likelihood function
    l = 0
    for i in (np.arange(n)+1):
        delta = (x[i-1] - xbar).reshape(pL,pR) - gamma2.dot(beta2).dot(fy[i-1]).dot(beta1.T).dot(gamma1.T)
        l += np.sum(temp.diagonal())

    L1 = l
    

# In[] 

# In[] 

# In[] 


# In[] 


# In[] return estimators
X11 = np.zeros((n,dL,dR))
for i in (np.arange(n)+1):
    X11[i-1] = gamma2.T.dot(np.mat(Mhat).I).dot(x[i-1]).dot(np.mat(omegahat).I).dot(gamma1)
est = [gamma1,gamma2,beta1,beta2,omegahat,Mhat,L1,X11]


# In[]


# In[]
# =============================================================================
# img1 = cv2.imread(filenames,cv2.IMREAD_GRAYSCALE)   #读取图片，第二个参数表示以灰度图像读入
# if img1 is None:                   #判断读入的img1是否为空，为空就继续下一轮循环
#             continue
# res1= cv2.resize(img1,(28,28))              #对图片进行缩放，第一个参数是读入的图片，第二个是制定的缩放大小
# res1_1 = res1.reshape(1,784)/255       #将表示图片的二维矩阵转换成一维
# res1_1_1 = res1_1.tolist()                     #将numpy.narray类型的矩阵转换成list
# train_set_x.append(res1_1_1)  
# =============================================================================
 
# In[]

# In[]


# In[]



# In[]


