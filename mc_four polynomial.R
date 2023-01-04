rm(list=ls())
set.seed(125)
T<-1000
n<-500
library(ggplot2)
library(dplyr)
library(rdd)
library(rddtools)
library(rddensity)
library(rdrobust)
##firstly,the data generating function:
DGF<-function(n,a){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x^2+x^4+10,x^2+x^4+1)
  w_s<-ifelse(x>=0,1,0)
  data<-data.frame(x,y,w_s)
  return(data)
}
mean_two<-matrix(NA,4,4)
##non-parametric method:
kernel <- c("triangular","epanechnikov","uniform")
## we will see the treatment effect in different kernel:
for (ker in 1:length(kernel)){
  a<-c()
  np<-matrix(NA,1000,4)
  for (i in 1:T){
    for (k in 1:4){
      a[i]<-runif(1,min=1,max=20)
      data<-DGF(n,a[i])
      locfit<-rdrobust(data$y,data$x,p=k,c=0,kernel=kernel[ker])
      np[i,k]<-locfit$coef[1]
    }
  }
  mean_two[ker,1:4] <- c(mean(np[1:1000,1]),mean(np[1:1000,2]),mean(np[1:1000,3]),mean(np[1:1000,4]))
}
mean_two
##for the parametric method:
p_b1<-c()
p_b2<-c()
p_b3<-c()
p_b4<-c()
c=0
for (i in 1:T){
  a[i]<-runif(1,min=1,max=50)
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
mean_two
#from this we can see that with the increasing of P, the estimation accuracy becomes higher