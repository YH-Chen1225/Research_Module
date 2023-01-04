set.seed(122)
library(ggplot2)
library(dplyr)
library(rdd)
library(rddtools)
library(rddensity)
library(rdrobust)

x<-runif(100,-2,2)
ran<- sample(0:100, 10, TRUE)
y<-ifelse(x>=0,3*x+10,3*x+1)#from here we know that the treatment affect should be 9

for (i in ran){
  if (x[i] <= 0){
    y[i] = 3*x[i]+10
  }
  else {y[i] = 3*x[i]+1}
}
plot(x, y, xlab = "x", ylab = "y", pch = 20, cex.axis = 1.5, cex.lab = 1.5)
abline(v = 0)

a <- rnorm(1,0,1)
a
#Data generating process With always taker and never taker
#Linear_FRD
DGF<-function(n,a){
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

#Non_Linear_FRD
DGF2 <- function(n,a){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x^2+50,-x^2+1)
  y<-y+rnorm(n,0,5)
  w_s<-ifelse(x>=0,1,0)
  ran<- sample(1:n, n/10, FALSE)
  for (i in ran){
    y[i] <- ifelse(x[i]>=0,x[i]^2+1+rnorm(1,0,5),-x[i]^2+50+rnorm(1,0,5))}
  data<-data.frame(x,y,w_s)
  return(data)
}

#Checking the data generating process
data <- DGF(100,10)
plot(data$x, data$y, xlab = "x", ylab = "y", pch = 20, cex.axis = 1.5, cex.lab = 1.5)
abline(v = 0)

data2 <- DGF2(100,10)
plot(data2$x, data2$y, xlab = "x", ylab = "y", pch = 20, cex.axis = 1.5, cex.lab = 1.5)
abline(v = 0)


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

#DGP can be DGF or DGF2, which is for non-linear simulation
m <- FRD(40,5,DGP=DGF,kernel=kernel,n=500)
m



