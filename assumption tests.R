set.seed(122133)
library(ggplot2)
library(dplyr)
library(rdd)
library(rddtools)
library(rddensity)
library(rdrobust)
########generate new data
DGF<-function(n,a){
  x<-runif(n,-a,a)
  y<-ifelse(x>=0,x^2+50,-x^2+1)
  y<-y+rnorm(n,0,5)
  w_s<-ifelse(x>=0,1,0)
  ran<- sort(sample(1:n, n/10, FALSE))
  for (i in ran){
    y[i] <- ifelse(x[i]>=0,x[i]^2+1+rnorm(1,0,5),-x[i]^2+50+rnorm(1,0,5))
    w_s[i]<-ifelse(x[i]>=0,0,1)
  }
  data<-data.frame(x,y,w_s)
  return(data)
}
data<-DGF(n=1000,a=20)
c<-0
x<-data$x
y<-data$y
w_s<-data$w_s
########1: show the Scatterplot:
plot(x, y, xlab = "x", ylab = "y", pch = 20,cex.axis = 1.5, cex.lab =1)
abline(v = 0)
########2:see the scatter is sharp or fuzzy
ggplot(data=data,mapping=aes(x=x,y=w_s,color=w_s))+
  geom_point(size=0.5,alpha=1,position = position_jitter(width=0,height=0.2,seed=1234))+
  geom_vline(xintercept = 0)
########3:check for the discontinuity in running variable around cutpoint:
data %>%
  group_by(w_s,x>=0) %>%
  summarise(count=n())
ggplot(data=data,mapping=aes(x=x))+###here miss a fill?
  geom_histogram(binwidth = 0.5,boundary=0,color="white")+
  labs(title = 'Histogram of the Rating varibale', y = 'Number of Observation', x = 'x') + 
  geom_vline(xintercept = 0)
#or there is another way to do this:
rdplotdensity(rdd=rddensity(data$x,c=0),X=data$x,type="both")
########4:placebo test
dat_rddtools<-rddtools::rdd_data(y=y,x=x,data=data,cutpoint=0)
llm_rddtools<-rddtools::rdd_reg_np(dat_rddtools)
fs<-rddtools::plotPlacebo(llm_rddtools)
summary(fs)
########5:check for discontinuity in outcome across running variable:
ggplot(data,aes(x=x,y=y,color=w_s))+
  geom_vline(xintercept = 0)+
  geom_point(size=0.75,alpha=0.5)+
  geom_smooth(data=filter(data,x<=0),method="loess")+
  geom_smooth(data=filter(data,x >0),method="loess")
#and we can see that there is a gap.so that is valid.








##########other notes,not important.
#Then:Fitting Graph and (Raw comparison of means):
rdplot(y, x, p = 1, col.lines = "red", col.dots = "black", title = "linear fitting",
       x.label = "x",y.label = "y", y.lim = c(-400,400))
rdplot(y, x, p = 2, col.lines = "red", col.dots = "black", title = "non-linear fitting",
       x.label = "x",y.label = "y", y.lim = c(-400,400))
#if we need to add the interval:
# rdplot(y, x, p = 2,ci=95, col.lines = "red", col.dots = "black", title = "non-linear fitting",
       # x.label = "x",y.label = "y", y.lim = c(-5,400))
#another way to create the sub-data:
#data2<-subset(data,data$x>-1&data$x< 1,select=x:y)
#local comparison of means:
# rdplot(y[abs(x) <= 1], x[abs(x) <= 1],p = 2,ci=95, col.lines = "red", col.dots = "black", title = "",
       # x.label = "x", y.label = "y", y.lim = c(-400,400))
##aggregate or \smooth" the data before plotting which is same with the first way:local regreesion of means.
#out = rdplot(y, x, nbins = c(20,20), binselect = 'esmv', y.lim = c(-5,15))
#summary(out)
########assumption 1: McCrary test
#from the graph, we can see that in cutoff, there is no obvious difference.
rdd::DCdensity(x,0,ext.out=T,plot=F)
#z=1.6,p=0.1,so no violation
########assumption 2 :placebo test
dat_rddtools<-rddtools::rdd_data(y=y,x=x,data=data,cutpoint=0)
llm_rddtools<-rddtools::rdd_reg_np(dat_rddtools)
#it is another way to show the assumption 1 use rddtools
#rddtools::dens_test(llm_rddtools)
fs<-rddtools::plotPlacebo(llm_rddtools)
summary(fs)

