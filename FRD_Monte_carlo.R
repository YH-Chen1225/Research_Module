set.seed(122)
library(ggplot2)
library(dplyr)
library(rdd)
library(rddtools)
library(rddensity)
library(rdrobust)
library(MuMIn)

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

#Data generating process With always taker and never taker
#Linear_FRD
DGF1<-function(n,a,gap){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,3*x+1+gap,3*x+1)
  y<-y+rnorm(n,0,1)
  w_s<-ifelse(x>=0,1,0)
  ran<- sort(sample(1:n, n/10, FALSE))
  for (i in ran){
    y[i] <- ifelse(x[i]>=0,3*x[i]+1+rnorm(1,0,1),3*x[i]+1+gap+rnorm(1,0,1))
    w_s[i]<-ifelse(x[i]>=0,0,1)
  }
  data<-data.frame(x,y,w_s)
  return(data)
}

#Non_Linear_FRD
DGF2 <- function(n,a,gap){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x^2+1+gap,-x^2+1)
  y<-y+rnorm(n,0,5)
  w_s<-ifelse(x>=0,1,0)
  ran<-sort(sample(1:n, n/10, FALSE))
  for (i in ran){
    y[i] <- ifelse(x[i]>=0,x[i]^2+1+rnorm(1,0,5),-x[i]^2+1+gap+rnorm(1,0,5))
    w_s[i]<-ifelse(x[i]>=0,0,1)
  }
  data<-data.frame(x,y,w_s)
  return(data)
}

#Third_order_FRD
DGF3 <- function(n,a,gap){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x^2+x^3+gap,x^2+x^3+1)
  y<-y+rnorm(n,0,20)
  w_s<-ifelse(x>=0,1,0)
  ran<- sort(sample(1:n, n/10, FALSE))
  for (i in ran){
    y[i] <- ifelse(x[i]>=0,x[i]^2+x[i]^3+1+rnorm(1,0,20),x[i]^2+x[i]^3+gap+rnorm(1,0,20))
    w_s[i]<-ifelse(x[i]>=0,0,1)
  }
  data<-data.frame(x,y,w_s)
  return(data)
}

#Fourth_order_FRD
DGF4<-function(n,a,gap){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x^2+x^3+x^4+1+gap,x^2+x^3+x^4+1)
  y<-y+rnorm(n,0,50)
  w_s<-ifelse(x>=0,1,0)
  ran<- sort(sample(1:n, n/10, FALSE))
  for (i in ran){
    y[i] <- ifelse(x[i]>=0,x[i]^2+x[i]^3+x[i]^4+1+rnorm(1,0,50),x[i]^2+x[i]^3+x[i]^4+1+gap+rnorm(1,0,50))
    w_s[i]<-ifelse(x[i]>=0,0,1)
  }
  data<-data.frame(x,y,w_s)
  return(data)
}


#Checking the data generating process
data <- DGF1(100,10,10)
plot(data$x, data$y, xlab = "x", ylab = "y", pch = 20, cex.axis = 1.5, cex.lab = 1.5)
abline(v = 0)

data2 <- DGF2(100,10,50)
plot(data2$x, data2$y, xlab = "x", ylab = "y", pch = 20, cex.axis = 1.5, cex.lab = 1.5)
abline(v = 0)

data3 <- DGF3(100,10,500)
plot(data3$x, data3$y, xlab = "x", ylab = "y", pch = 20, cex.axis = 1.5, cex.lab = 1.5)
abline(v = 0)


#Preparation for Monte Carlo
kernel <- c("triangular","epanechnikov","uniform")

#Monte Carlo simulation
FRD <- function(T,order,n,DGP,kernel,gap){
  mean_kernel <- matrix(NA,length(kernel),order)  
  for (ker in 1:length(kernel)){
    a<-c()
    np<-matrix(NA,T,order)
    for (i in 1:T){
      for (k in 1:order){
        a[i]<-runif(1,min=1,max=10)
        data<-DGP(n,a[i],gap)
        locfit<-rdrobust(data$y,data$x,fuzzy=data$w_s,p=k,c=0,kernel=kernel[ker])
        np[i,k]<-locfit$coef[1]
      }
  }
    for (z in 1:order){
    mean_kernel[ker,z] <- mean((np[1:T,z]-gap)^2)}
  }
return(mean_kernel)
}

#DGP can be DGF ,DGF2 and DGF3, which is for non-linear simulation
m_1_1 <- FRD(1000,4,DGP=DGF1,kernel=kernel,n=700,gap = 1)#n cannot be too small
m_1_1
mean_1_1 <- apply(m_1_1,2,mean)
mean_1_1

m_1_2 <- FRD(1000,4,DGP=DGF1,kernel=kernel,n=800,gap = 0.1)#n cannot be too small
m_1_2
mean_1_2 <- apply(m_1_2,2,mean)
mean_1_2

m_1_3 <- FRD(1000,4,DGP=DGF1,kernel=kernel,n=800,gap = 0.01)#n cannot be too small
m_1_3
mean_1_3 <- apply(m_1_3,2,mean)
mean_1_3


m_1_4 <- FRD(1000,4,DGP=DGF1,kernel=kernel,n=800,gap = 0.001)#n cannot be too small
m_1_4
mean_1_4 <- apply(m_1_4,2,mean)
mean_1_4

m_1_5 <- FRD(1000,4,DGP=DGF1,kernel=kernel,n=800,gap = 0.0001)#n cannot be too small
m_1_5
mean_1_5 <- apply(m_1_5,2,mean)
mean_1_5

m_1_6 <- FRD(1000,4,DGP=DGF1,kernel=kernel,n=800,gap = 0.00001)#n cannot be too small
m_1_6
mean_1_6 <- apply(m_1_6,2,mean)
mean_1_6


###Parametric Method
vars <- c("w_s","data$x_centered","I(x_centered^2)","I(x_centered^3)","I(x_centered^4)")
para <- function(T,vars,c,DGP,gap,n){
  a <- c()
  coe <- matrix(NA,T,4)
  for (i in 1:T){
    a[i]<-runif(1,min=1,max=10)
    data<-DGP(n,a[i],gap)
    data$x_centered<-data$x-c
    for (z in 1:(length(vars)-1)){
      formula <- as.formula(paste("y",paste(vars[0:z+1],collapse = "+"),sep = "~"))
      fit <- lm(formula,data=data)
      coe[i,z] <- fit$coefficients[2]
    }
  }
  sqr <- (coe-gap)^2
  mse <- apply(sqr,2,mean)
  return(mse)
}


apply(coe, 2, mean)
##DGP can be DGF or DGF2 or DGF3, which is for non-linear data generating process
mse <- para(1000,vars=vars,c=0,DGP=DGF1,gap = 1,n=800)
mse

mse2 <- para(1000,vars=vars,c=0,DGP=DGF1,gap = 0.1,n=800)
mse2

mse3 <- para(1000,vars=vars,c=0,DGP=DGF1,gap = 0.01,n=800)
mse3

mse4 <- para(1000,vars=vars,c=0,DGP=DGF1,gap = 0.001,n=800)
mse4

mse5 <- para(1000,vars=vars,c=0,DGP=DGF1,gap = 0.0001,n=800)
mse5


#Simple testing
data <- DGF3(100,30)
fit <- lm(y~w_s+x,data=data)
fit$coefficients

ma <- matrix(4,2,2)
q <- ma-1
q_2 <- q^2
q_2

