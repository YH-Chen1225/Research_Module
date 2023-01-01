set.seed(122)
library(ggplot2)
library(dplyr)
install.packages('rdd')
install.packages('rddtools')
library(rdd)
library(rddtools)
install.packages('rddensity')
library(rddensity)
install.packages('rdrobust')
library(rdrobust)
########generate new data
x<-runif(100,-2,2)
y<-ifelse(x>=0,3*x+10,2*x+1)
w_s<-ifelse(x>=0,1,0)
data<-data.frame(x,y,w_s)
#Covariance of the forcing variable and outcome
#sigma12 <- corr * sqrt(sigma1) * sqrt(sigma2)
# Mean of prognosis and outcome
#mu<- c(10,90)
# Covariance Matrix
#sigma <- matrix(c(sigma1, sigma12, sigma12, sigma2), nrow = 2, byrow = TRUE)
#dataset <- rmvnorm(100, mean = mu, sigma = sigma)
#dataset        <- as.data.frame(dataset)
#names(dataset) <- c("forcing", "outcome")
########set the cut-point
c<-0

########Firstly, the Scatterplot:
plot(x, y, xlab = "x", ylab = "y", pch = 20, cex.axis = 1.5, cex.lab = 1.5)
abline(v = 0)
########then see the scatter is sharp or fuzzy
ggplot(data=data,mapping=aes(x=x,y=w_s,color=w_s))+
  geom_point(size=0.5,alpha=1,position = position_jitter(width=0,height=0.1,seed=1234))+
  geom_vline(xintercept = 0)
# Then:Fitting Graph and (Raw comparison of means):
rdplot(y, x, p = 1, col.lines = "red", col.dots = "black", title = "linear fitting",
       x.label = "x",y.label = "y", y.lim = c(-5,20))
rdplot(y, x, p = 2, col.lines = "red", col.dots = "black", title = "non-linear fitting",
       x.label = "x",y.label = "y", y.lim = c(-5,20))
#if we need to add the interval:
rdplot(y, x, p = 2,ci=95, col.lines = "red", col.dots = "black", title = "non-linear fitting",
       x.label = "x",y.label = "y", y.lim = c(-5,20))
#another way to create the sub-data:
#data2<-subset(data,data$x>-1&data$x< 1,select=x:y)
#local comparison of means:
rdplot(y[abs(x) <= 1], x[abs(x) <= 1], nbins = c(2500, 500), p = 4, col.lines = "red", col.dots = "black", title = "",
       x.label = "x", y.label = "y", y.lim = c(-5,20))
##aggregate or \smooth" the data before plotting which is same with the first way:local regreesion of means.
#out = rdplot(y, x, nbins = c(20,20), binselect = 'esmv', y.lim = c(-5,15))
#summary(out)


########assumption 1: McCrary test
#Firstly, we can see the plot:
ggplot(data=data,aes(x=x))+###here miss a fill?
  geom_histogram(binwidth = 0.1,boundary=0,color="white")+
  geom_vline(xintercept = 0)
#from the graph, we can see that in cutoff, there is no obvious difference.
rdd::DCdensity(x,0,ext.out=T,plot=F)
#z=1.6,p=0.1,so no violation
#there is another way:
rdplotdensity(rdd=rddensity(data$x,c=0),X=data$x,type="both")
##finding discontinuity in outcome##maybe the assumption 2?
ggplot(data,aes(x=x,y=y,color=w_s))+
  geom_vline(xintercept = 0)+
  geom_point(size=0.75,alpha=0.5)+
  geom_smooth(data=filter(data,x<=0),method="lm")+
  geom_smooth(data=filter(data,x >0),method="lm")
########assumption 2 :placebo test
dat_rddtools<-rddtools::rdd_data(y=y,x=x,data=data,cutpoint=0)
llm_rddtools<-rddtools::rdd_reg_np(dat_rddtools)
#it is another way to show the assumption 1 use rddtools
#rddtools::dens_test(llm_rddtools)
fs<-rddtools::plotPlacebo(llm_rddtools)
summary(fs)

########estimate the effect:
summary(rdbwselect(y=y,x=x,all=TRUE))
#we also can use robust packet to select automatically:
locfit<-rdrobust(y,x,c=0,p=1,kernel='triangular',bwselect='mserd')
summary(locfit)
#to plot the rdrobust
rdplot(y,x,c=0,p=1,kernel='triangular')
rdplot(y,x,c=0,p=1,kernel = 'epanechnikov')
#from the result,we can see that the best bandwidth is 1.037,here the kernel
#is triangular,we can choose another one.
#####non-parametric Local Polynomial Regression
locfit2<-rdrobust(y,x,c=0,p=2,kernel='triangular',bwselect='mserd')
summary(locfit2)
locfit21<-rdrobust(y,x,c=0,p=2,kernel='epanechnikov',bwselect='mserd')
summary(locfit21)
locfit3<-rdrobust(y,x,c=0,p=3,kernel='triangular',bwselect='mserd')
summary(locfit3)
#####Parametric polynomial regression:different situations:
data$x_del<-data$x-c
fit1<-lm(y~x_del+w_s,data=data)
summary(fit1)
fit2<-lm(y~x_del+I(x_del^2)+w_s,data=data)
summary(fit2)
fit3<-lm(y~x_del+I(x_del^2)+I(x_del^3)+w_s,data=data)
summary(fit3)
##narrow the range  of x :
fit11<-lm(y~x_del+w_s,data=filter(data,x>=-1&x<=1))
summary(fit11)
