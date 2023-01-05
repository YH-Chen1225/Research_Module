set.seed(122)
library(ggplot2)
library(dplyr)
library(rdd)
library(rddtools)
library(rddensity)
library(rdrobust)
library(MuMIn)
DGF1<-function(n,a){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,3*x+10,3*x+1)
  y<-y+rnorm(n,0,1)
  w_s<-ifelse(x>=0,1,0)
  ran<- sample(1:n, n/10, FALSE)
  for (i in ran){
    y[i] <- ifelse(x[i]>=0,3*x[i]+1+rnorm(1,0,1),3*x[i]+10+rnorm(1,0,1))}
  data<-data.frame(x,y,w_s)
  return(data)
}

#Non_Linear_FRD(Gap should be 9)
DGF2 <- function(n,a){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x^2+10,-x^2+1)
  y<-y+rnorm(n,0,5)
  w_s<-ifelse(x>=0,1,0)
  ran<- sample(1:n, n/10, FALSE)
  for (i in ran){
    y[i] <- ifelse(x[i]>=0,x[i]^2+1+rnorm(1,0,5),-x[i]^2+10+rnorm(1,0,5))}
  data<-data.frame(x,y,w_s)
  return(data)
}

#Third_order_FRD(Gap should be 9)
DGF3 <- function(n,a){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x^2+x^3+10,x^2+x^3+1)
  y<-y+rnorm(n,0,20)
  w_s<-ifelse(x>=0,1,0)
  ran<- sample(1:n, n/10, FALSE)
  for (i in ran){
    y[i] <- ifelse(x[i]>=0,x[i]^2+x[i]^3+1+rnorm(1,0,20),x[i]^2+x[i]^3+10+rnorm(1,0,20))}
  data<-data.frame(x,y,w_s)
  return(data)
}
DGF4<-function(n,a){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x^2+x^3+x^4+10,x^2+x^3+x^4+1)
  y<-y+rnorm(n,0,50)
  w_s<-ifelse(x>=0,1,0)
  ran<- sample(1:n, n/10, FALSE)
  for (i in ran){
    y[i] <- ifelse(x[i]>=0,x[i]^2+x[i]^3+x[i]^4+1+rnorm(1,0,50),x[i]^2+x[i]^3+x[i]^4+10+rnorm(1,0,50))}
  data<-data.frame(x,y,w_s)
  return(data)
}
#Preparation for Monte Carlo
kernel <- c("triangular","epanechnikov","uniform")

#Monte Carlo simulation
FRD <- function(T,order,n,DGP,kernel){
  mean_kernel <- matrix(NA,length(kernel),order)  
  for (ker in 1:length(kernel)){
    a<-c()
    np<-matrix(NA,T,order)
    for (i in 1:T){
      for (k in 1:order){
        a[i]<-runif(1,min=1,max=10)
        data<-DGP(n,a[i])
        locfit<-rdrobust(data$y,data$x,fuzzy=data$w_s,p=k,c=0,kernel=kernel[ker])
        np[i,k]<-locfit$coef[1]
      }
    }
    for (z in 1:order){
      mean_kernel[ker,z] <- mean(np[1:T,z])}
  }
  return(mean_kernel)
}

#DGP can be DGF1,DGF2,DGF3 and DGF4.
#then we can use the above fuction to calculate the treatment effect in different kernel:
T<-1000
p<-4
n<-300
#I hope I can think a way to make the below codes simpler.
mean1<-matrix(NA,3,p)#mean1 is the mean value of the treatment effect in 1-polynomial 
for (i in 1:3){
  mean1[i,]<-FRD(T,p,n,DGF1,kernel[i])
}
mean2<-matrix(NA,3,p)#mean2 is the mean value of the treatment effect in 2-polynomial 
for (i in 1:3){
  mean2[i,]<-FRD(T,p,n,DGF2,kernel[i])
}
mean3<-matrix(NA,3,p)#mean3 is the mean value of the treatment effect in 3-polynomial 
for (i in 1:3){
  mean3[i,]<-FRD(T,p,n,DGF3,kernel[i])
}
mean4<-matrix(NA,3,p)#mean4 is the mean value of the treatment effect in 1-polynomial 
for (i in 1:3){
  mean4[i,]<-FRD(T,p,n,DGF4,kernel[i])
}
###Parametric Method
vars <- c("w_s","data$x_centered","I(x_centered^2)","I(x_centered^3)","I(x_centered^4)")

para <- function(T,vars,c,DGP){
  a<-c()
  coe <- matrix(NA,T,4)
  for (i in 1:T){
    a[i]<-runif(1,min=1,max=10)
    data<-DGP(n,a[i])
    data$x_centered<-data$x-c
    for (z in 1:(length(vars)-1)){
      formula <- as.formula(paste("y",paste(vars[0:z+1],collapse = "+"),sep = "~"))
      fit <- lm(formula,data=data)
      coe[i,z] <- fit$coefficients[2]
    }
  }
  avg <- apply(coe, 2, mean)
  return(avg)
}
avg<-para(T,vars=vars,c=0,DGP=DGF1)
avg
