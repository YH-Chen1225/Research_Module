rm(list=ls())
set.seed(122122)
library(ggplot2)
library(dplyr)
library(rdd)
library(rddtools)
library(rddensity)
library(rdrobust)
library(MuMIn)
DGF1<-function(n,a,gap){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x+1+gap,x+1)
  y<-y+rnorm(n,0,1)
  w_s<-ifelse(x>=0,1,0)
  data<-data.frame(x,y,w_s)
  return(data)
}
# data<-DGF1(100,10,9)
# plot(data$x, data$y, xlab = "x", ylab = "y", pch = 20, cex.axis = 1.5, cex.lab = 1.5,
#      main="DGF1(n=100,a=10,gap=9)"
#      )
# abline(v = 0)
DGF2 <- function(n,a,gap){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x^2+1+gap,-x^2+1)
  y<-y+rnorm(n,0,5)
  w_s<-ifelse(x>=0,1,0)
  data<-data.frame(x,y,w_s)
  return(data)
}
# data2<-DGF2(100,10,99)
# plot(data2$x, data2$y, xlab = "x", ylab = "y", pch = 20, cex.axis = 1.5, cex.lab = 1.5,
# main="DGF2(n=100,a=10,gap=99)")
# abline(v = 0)

DGF3<-function(n,a,gap){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x^2+x^3+1+gap,x^2+x^3+1)
  y<-y+rnorm(n,0,20)
  w_s<-ifelse(x>=0,1,0)
  data<-data.frame(x,y,w_s)
  return(data)
}
# data3<-DGF3(100,10,199)
# plot(data3$x, data3$y, xlab = "x", ylab = "y", pch = 20, cex.axis = 1.5, cex.lab = 1.5,
#      main="DGF3(n=100,a=10,gap=199)")
# abline(v = 0)
DGF4<-function(n,a,gap){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x^2+x^3+x^4+1+gap,-x^2+x^3-x^4+1)
  y<-y+rnorm(n,0,50)
  w_s<-ifelse(x>=0,1,0)
  data<-data.frame(x,y,w_s)
  return(data)
}
 # data4<-DGF4(100,10,999)
 # plot(data4$x, data4$y, xlab = "x", ylab = "y", pch = 20, cex.axis = 1.5, cex.lab = 1.5,
 # main="DGF4(n=100,a=10,gap=999)")
 # abline(v = 0)

##non-parametric method:
kernel <- c("triangular","epanechnikov","uniform")
## we will see the treatment effect in different kernel:
SRD <- function(T,order,n,DGP,kernel,gap){
  mse_kernel <- matrix(NA,length(kernel),order)  
  for (ker in 1:length(kernel)){
    a<-c()
    np<-matrix(NA,T,order)
    for (i in 1:T){
      for (k in 1:order){
        a[i]<-20
        data<-DGP(n,a[i],gap)
        locfit<-rdrobust(data$y,data$x,p=k,c=0,kernel=kernel[ker])
        np[i,k]<-locfit$coef[1]
      }
    }
    for (z in 1:order){
      mse_kernel[ker,z] <- mean((np[1:T,z]-gap)^2)#
                          #sqrt(mean((np[1:T,z]-gap)^2))/mean(np[1:T,z])#NRMSE
      }
  }
  return(mse_kernel)
}
m_s1_1 <- SRD(1000,4,DGP=DGF1,kernel=kernel,n=1000,gap = 9)#n cannot be too small
m_s1_1

m_s1_2 <- SRD(1000,4,DGP=DGF2,kernel=kernel,n=1000,gap = 9)#n cannot be too small
  m_s1_2

m_s1_3 <- SRD(1000,4,DGP=DGF3,kernel=kernel,n=1000,gap = 9)#n cannot be too small
m_s1_3

m_s1_4 <- SRD(1000,4,DGP=DGF4,kernel=kernel,n=1000,gap = 9)#n cannot be too small
m_s1_4
###we want to check that when n is bigger, if the local linear estimation is better:
n_s<-seq(1000,10000,by=1000)
p=4
# mse_ns_t<-matrix(NA,length(n_s),4)
mse_ns_e<-matrix(NA,length(n_s),4)
# mse_ns_u<-matrix(NA,length(n_s),4)
mse_ns_e[10,]<-SRD(1000,p,10000,DGF4,kernel[2],gap=9)
for (i in 1:length(n_s)){
  # mse_ns_t[i,]<-SRD(1000,p,n_s[i],DGF1,kernel[1],gap=9)
  mse_ns_e[i,]<-SRD(1000,p,n_s[i],DGF4,kernel[2],gap=9)
  # mse_ns_u[i,]<-SRD(1000,p,n_s[i],DGF1,kernel[3],gap=9)
}
########find the critical value of the gap:
gap_s<-c(1,0.1,0.01,0.001,0.0001,0.00001,0.000001,0.0000001,0.00000001,0.000000001,0.0000000001)
##for DGF1:
m_s2_1<-matrix(NA,length(gap_s),4)
m_s2_2s<-matrix(NA,length(gap_s),4)
m_s2_3<-matrix(NA,length(gap_s),4)
m_s2_4<-matrix(NA,length(gap_s),4)
for (i in 1:length(gap_s)){
  m_s2_2s[i,]=SRD(1000,4,DGP=DGF2,kernel=kernel[2],n=1000,gap = gap_s[i])#n cannot be too small
} 


########parametric method:
vars <- c("w_s","data$x_centered","I(data$x_centered^2)","I(data$x_centered^3)","I(data$x_centered^4)")
para <- function(T,vars,c,DGP,gap,n){
  a<-c()
  coe <- matrix(NA,T,4)
  for (i in 1:T){
    a[i]<-runif(1,min=1,max=20)
    data<-DGP(n,a[i],gap)
    data$x_centered<-data$x-c
    for (z in 1:(length(vars)-1)){
      formula <- as.formula(paste("y",paste(vars[0:z+1],collapse = "+"),sep = "~"))
      fit <- lm(formula,data=data)
      coe[i,z] <- fit$coefficients[2]
    }
  }
  sqr <- (coe-gap)^2
  mse<-apply(sqr,2,mean)
  return(mse)
}
avg<-matrix(NA,4,4) 
avg[1,]<- para(T,vars=vars,c=0,DGP=DGF1,gap=9,n=1000)
avg[2,]<- para(T,vars=vars,c=0,DGP=DGF2,gap=9,n=1000)
avg[3,]<- para(T,vars=vars,c=0,DGP=DGF3,gap=9,n=1000)
avg[4,]<- para(T,vars=vars,c=0,DGP=DGF4,gap=9,n=1000)


