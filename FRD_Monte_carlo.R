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
##########################FRD_Linear
#With always taker and never taker
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

#Checking the data generating process
data <- DGF(100,10)
plot(data$x, data$y, xlab = "x", ylab = "y", pch = 20, cex.axis = 1.5, cex.lab = 1.5)
abline(v = 0)

#Monte Carlo Simulation
T<-1000
a<-c()
n<-500
mean_kernel <- matrix(NA,3,4)
kernel <- c("triangular","epanechnikov","uniform")

for (ker in 1:length(kernel)){
  a<-c()
  np<-matrix(NA,1000,4)
  for (i in 1:T){
    for (k in 1:4){
    a[i]<-runif(1,min=1,max=10)
    data<-DGF(n,a[i])
    locfit<-rdrobust(data$y,data$x,fuzzy=data$w_s,p=k,c=0,kernel=kernel[ker])
    np[i,k]<-locfit$coef[1]
    }
  }
  mean_kernel[ker,1:4] <- c(mean(np[1:1000,1]),mean(np[1:1000,2]),mean(np[1:1000,3]),mean(np[1:1000,4]))
}

mean_kernel




