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

decomp(a,-1)
# In[]

# In[]  read images
n=40
array_of_img = []
for i in ["glass/glass","paper/paper","metal/metal"]:
    for j in (np.arange(n)+1):
        img = mpimg.imread(i + str(j) + ".jpg")
        img_grey = np.dot(img[...,:3], [0.299, 0.587, 0.114])  # convert rgb to grey
        array_of_img.append(img_grey)

array_of_img[0].shape
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
for i in (np.arange(n)+1):
    covx += array_of_img[i].T.dot(array_of_img[i])

gamma1 = np.linalg.eig(covx)[1][:,:dR]  # generate gamma1
beta1 = np.array([r.random() for i in np.arange(dR*rR)]).reshape((dR,rR))  # generate beta1
omeghat = np.identity(pR)*[np.abs(r.random()) for i in np.arange(pR)]  # generate omegahat
Mhat = np.identity(pL)*[np.abs(r.random()) for i in np.arange(pL)]  # generate Mhat
# In[]
x = array_of_img
# standardizes the predictors
z = [decomp(Mhat,-1/2).dot(x[i]).dot(decomp(omeghat,-1/2)) for i in (np.arange(n)+1)]

# In[]  update gamma2 and beta2
gamma1z = decomp(omeghat,-1/2).dot(gamma1)
# In[]
XL = []
for i in (np.arange(n)):
    XL[((i-1)*dR+1) : (i*dR+1)] = (z[i].dot(gamma1z)).T

#  XL[((i-1)*dR+1) : (i*dR),] <- t(z[,,i] %*% gamma1z)
#      FL[((i-1)*dR+1) : (i*dR),] <- t(fy[,,i] %*% t(beta1))
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
 
# In[] standardizes the predictors

np.dot(decomp(Mhat,-1/2), x[,,i]).dot(decomp(omegahat,-1/2))



# In[]

mat = np.matrix([[1,2],[3,4]])

np.dot(mat_vec.T, mat_vec)
#np.dot(mat_vec,)
np.diagonal(mat)
a=np.matrix([[1,2],[3,4]])

# In[]
mat_val
# In[]

# In[]


# In[]



# In[]


