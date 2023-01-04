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

##########################FRD_Linear
#With always taker and non takere
DGF<-function(n,a){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,3*x+10,3*x+1)
  w_s<-ifelse(x>=0,1,0)
  ran<- sample(1:n-1, n/10, FALSE)
  for (i in ran){
    y[i] <- ifelse(x[i]>=0,3*x[i]+1,3*x[i]+10)}
  data<-data.frame(x,y,w_s)
  return(data)
}

#Only compiler included
DGF<-function(n,a){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,3*x+10,2*x+1)
  w_s<-ifelse(x>=0,1,0)
  data<-data.frame(x,y,w_s)
  return(data)
}
#
######
data <- DGF(100,2)
locfit<-rdrobust(data$y,data$x,fuzzy=data$w_s,p=k,c=0,kernel = "uniform")
locfit$coef[1]




T<-1000
a<-c()
mean_kernel <- matrix(NA,3,4)
kernel <- c("triangular","epanechnikov","uniform")

for (ker in 1:length(kernel)){
  a<-c()
  np<-matrix(NA,1000,4)
  for (i in 1:T){
    for (k in 1:4){
    a[i]<-runif(1,min=1,max=20)
    data<-DGF(100,a[i])
    locfit<-rdrobust(data$y,data$x,fuzzy=data$w_s,p=k,c=0,kernel=kernel[ker])
    np[i,k]<-locfit$coef[1]
    }
  }
  mean_kernel[ker,1:4] <- c(mean(np[1:1000,1]),mean(np[1:1000,2]),mean(np[1:1000,3]),mean(np[1:1000,4]))
}

mean_kernel




