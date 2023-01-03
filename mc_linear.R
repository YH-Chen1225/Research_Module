rm(list=ls())
set.seed(122)
T<-1000
n<-300
library(ggplot2)
library(dplyr)
library(rdd)
library(rddtools)
library(rddensity)
library(rdrobust)
##firstly,the data generating function:
DGF<-function(n,a){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,3*x+10,2*x+1)
  w_s<-ifelse(x>=0,1,0)
  data<-data.frame(x,y,w_s)
  return(data)
}
mean_two<-matrix(NA,4,4)
##non-parametric method:
## we will see the treatment effect when the kernel is triangular:
a<-c()
b1_t<-c()
b2_t<-c()
b3_t<-c()
b4_t<-c()
for (i in 1:T){
  a[i]<-runif(1,min=1,max=20)
  data<-DGF(n,a[i])
  locfit1<-rdrobust(data$y,data$x,p=1,c=0,kernel="triangular")
  b1_t[i]<-locfit1$coef[1]
  locfit2<-rdrobust(data$y,data$x,p=2,c=0,kernel="triangular")
  b2_t[i]<-locfit2$coef[1]
  locfit3<-rdrobust(data$y,data$x,p=3,c=0,kernel="triangular")
  b3_t[i]<-locfit3$coef[1]
  locfit4<-rdrobust(data$y,data$x,p=4,c=0,kernel="triangular")
  b4_t[i]<-locfit4$coef[1]
}
b1_t
b2_t
b3_t
b4_t
mean_two[1,]<-c(mean(b1_t),mean(b2_t),mean(b3_t),mean(b4_t))
## we will see the treatment effect when the kernel is epanechnikov:
b1_e<-c()
b2_e<-c()
b3_e<-c()
b4_e<-c()
for (i in 1:T){
  a[i]<-runif(1,min=1,max=20)
  data<-DGF(n,a[i])
  locfit1<-rdrobust(data$y,data$x,p=1,c=0,kernel="epanechnikov")
  b1_e[i]<-locfit1$coef[1]
  locfit2<-rdrobust(data$y,data$x,p=2,c=0,kernel="epanechnikov")
  b2_e[i]<-locfit2$coef[1]
  locfit3<-rdrobust(data$y,data$x,p=3,c=0,kernel="epanechnikov")
  b3_e[i]<-locfit3$coef[1]
  locfit4<-rdrobust(data$y,data$x,p=4,c=0,kernel="epanechnikov")
  b4_e[i]<-locfit4$coef[1]
}
b1_e
b2_e
b3_e
b4_e
mean_two[2,]<-c(mean(b1_e),mean(b2_e),mean(b3_e),mean(b4_e))
## we will see the treatment effect when the kernel is uniform:
b1_u<-c()
b2_u<-c()
b3_u<-c()
b4_u<-c()
for (i in 1:T){
  a[i]<-runif(1,min=1,max=20)
  data<-DGF(n,a[i])
  locfit1<-rdrobust(data$y,data$x,p=1,c=0,kernel="uniform")
  b1_u[i]<-locfit1$coef[1]
  locfit2<-rdrobust(data$y,data$x,p=2,c=0,kernel="uniform")
  b2_u[i]<-locfit2$coef[1]
  locfit3<-rdrobust(data$y,data$x,p=3,c=0,kernel="uniform")
  b3_u[i]<-locfit3$coef[1]
  locfit4<-rdrobust(data$y,data$x,p=4,c=0,kernel="uniform")
  b4_u[i]<-locfit4$coef[1]
}
b1_u
b2_u
b3_u
b4_u
mean_two[3,]<-c(mean(b1_u),mean(b2_u),mean(b3_u),mean(b4_u))
##for the parametric method:
p_b1<-c()
p_b2<-c()
p_b3<-c()
p_b4<-c()
c=0
for (i in 1:T){
  a[i]<-runif(1,min=1,max=20)
  data<-DGF(n,a[i])
  data$x_centered<-data$x-c
  fit1<-lm(data$y~data$x_centered+data$w_s,data=data)
  p_b1[i]<-fit1$coefficients[3]
  fit2<-lm(data$y~data$x_centered+I(x_centered^2)+w_s,data=data)
  p_b2[i]<-fit2$coefficients[4]
  fit3<-lm(data$y~data$x_centered+I(x_centered^2)+I(x_centered^3)+w_s,data=data)
  p_b3[i]<-fit3$coefficients[5]
  fit4<-lm(data$y~data$x_centered+I(x_centered^2)+I(x_centered^3)+I(x_centered^4)+w_s,data=data)
  p_b4[i]<-fit4$coefficients[6]
}
p_b1
p_b2
p_b3
p_b4
mean_two[4,]<-c(mean(p_b1),mean(p_b2),mean(p_b3),mean(p_b4))
show(mean_two)
#from this we can see that with the increasing of P, the estimation accuracy becomes higher