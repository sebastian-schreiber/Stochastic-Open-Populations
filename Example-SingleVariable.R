# The Example-SingleVariable.R file for the Methods in Ecology and Evolution (2018) paper 
# "The structured demography of open populations in fluctuating environments" 
# by Sebastian Schreiber and Jacob Moore. 
# See the READ-ME.txt for an overview of the files and their dependencies. 

# This file studies the simple, scalar model 
# n(t+1)=exp(mua+sigmaa*Z1(t+1))*n(t)+exp(mub+sigmab*Z2(t+1))
# where Z1, Z2 are standard normals. 



# load packages
library(mvtnorm) # for the multivariate normal command
library(viridis)

# numerical simulator for the process
sim=function(parms,Tmax=100){
  set.seed(1)
  mua=parms$mua
  mub=parms$mub
  sigmaa=parms$sigmaa
  sigmab=parms$sigmab
  rho=parms$rho
  n=numeric(Tmax)
  n[1]=exp(mub+sigmab^2/2)/(1-exp(mua+sigmaa^2/2))
  Sigma=rbind(c(sigmaa^2,rho*sigmaa*sigmab),c(rho*sigmaa*sigmab,sigmab^2))
  for(t in 2:Tmax){
    X=exp(rmvnorm(n=1,mean=c(mua,mub),sigma=Sigma))
    n[t]=X[1]*n[t-1]+X[2]
  }
  return(n)
}


# analytic computation of variance for stationary distribution
v.func=function(parms){
  mua=parms$mua
  mub=parms$mub
  sigmaa=parms$sigmaa
  sigmab=parms$sigmab
  rho=parms$rho
  nbar=exp(mub+sigmab^2/2)/(1-exp(mua+sigmaa^2/2))
  abar=exp(mua+sigmaa^2/2)
  a2bar=exp(2*mua+2*sigmaa^2)
  bbar=exp(mub+sigmab^2/2)
  var=nbar^2*abar^2*(exp(sigmaa^2)-1)+bbar^2*(exp(sigmab^2)-1)+
    2*nbar*abar*bbar*(exp(sigmaa*sigmab*rho)-1)
  var=var/(1-a2bar)
  return(var)
}

# Figure code
pdf("single-variableA.pdf",width=7,height=4)
par(cex.lab=1.25,cex.axis=1.25,mar=c(4.1,4.1,1,0))
layout(rbind(c(1,2,3),c(4,4,4)))
hist.col='white'

# base parameters
parms=list(mua=-0.05,sigmaa=0.05,mub=0.1,sigmab=0.05,rho=0)
#critical sig values for finite mean and finite variance
critical.sig=sqrt(2*(-parms$mua))
critical.sig2=sqrt(-parms$mua)

#Histogram plots
# length of run for each histogram
Tmax=10000

# histogram 1
parms$sigmaa=0.01
out=sim(parms,Tmax)
print(var(out))
print(v.func(parms))
breaks=seq(0,max(out),length=40)
hist(out,breaks=breaks,freq=FALSE,main="",col=hist.col, xlab="population density",ylab="frequency")
legend("topleft",c(expression(sigma[a]==0.01),expression(bar(a)==0.951),expression(mu[a]-sigma[a]^2/2==0.05)),bty="n")

#histogram 2
parms$sigmaa=0.1
out=sim(parms,Tmax)
print(var(out))
print(v.func(parms))
breaks=seq(0,max(out),length=40)
hist(out,breaks=breaks,freq=FALSE,main="",col=hist.col, xlab="population density",ylab="")
legend("topright",c(expression(sigma[a]==0.1),expression(bar(a)==0.955),expression(mu[a]-sigma[a]^2/2==0.045)),bty="n")

#histogram 3
parms$sigmaa=0.2
print(exp(parms$mua+parms$sigmaa^2/2))
out=sim(parms,Tmax)
print(var(out))
print(v.func(parms))
breaks=seq(0,max(out),length=40)
hist(out,breaks=breaks,freq=FALSE,main="",col=hist.col, xlab="population density",ylab="")
legend("topright",c(expression(sigma[a]==0.2),expression(bar(a)==0.970),expression(mu[a]-sigma[a]^2/2==0.03)),bty="n")


# ploting means and variances as a function of sigmaa

# function to compute the mean
M=function(sigmaa)exp(parms$mub+parms$sigmab^2/2)/(1-exp(parms$mua+sigmaa^2/2))
# function to compute the SD
SD=function(sigmaa){
  parms$sigmaa=sigmaa
  return(sqrt(v.func(parms)))
}
# function to compute mean and SD numerically
numericMSD=function(x,Tmax){
  parms$sigmaa=x
  out=sim(parms,Tmax)
  return(c(mean(out),sd(out)))
}

# range of sigma values for analytic and numeric, respectively
sigs=seq(0,critical.sig,length=100)
sigs2=seq(0,critical.sig,by=0.05)
L=length(sigs2)

# compute the numerics
temp=function(x)numericMSD(x,10000)
numerics=sapply(sigs2,temp)
# improve numeric for 0.2,0.25,0.3	
longer=numericMSD(0.2,1000000)
numerics[,5]=longer
longer=numericMSD(0.25,1000000)
numerics[,6]=longer
longer=numericMSD(0.3,1000000)
numerics[,7]=longer

plot(sigs,M(sigs),lwd=2,type="l",xlab=expression(sigma[a]),
     ylab="mean and SD",ylim=c(0,200))
points(sigs2,numerics[1,],pch=4)
abline(v=critical.sig,lty=2)
sigs4=seq(0,critical.sig2*0.99999999,length=100)
second.col = '#969696'
lines(sigs4,sapply(sigs4,SD),lwd=2,col= second.col)
points(sigs2[1:5],numerics[2,1:5],pch=4,col= second.col)
abline(v=critical.sig2,lty=2,col= second.col)
legend('topleft', legend=c('mean', 'SD'), col=c('black', second.col), lty=c(1,1))
dev.off()





# # # plots of how the variance varies with rho (no effect on mean)
# # and variance

# parms=list(mua=-0.05,sigmaa=0.1,mub=0.1,sigmab=0.05,rho=0)



# k=25
# Tmax=1000
# sigs=seq(0,0.2,length=k)
# rhos=seq(0,2,length=k)
# std=matrix(0,k,k)
# c=std
# for(i in 1:k){
  # for(j in 1:k){
    # parms$sigmaa=sigs[i]
    # parms$sigmab=rhos[j]
    # std[i,j]=sqrt(v.func(parms))
    # out=sim(parms,Tmax)
    # out=closed(parms,out)
    # c[i,j]=mean(out)
  # }
# }

# pdf("single-variableB.pdf",width=6,height=5)
# par(cex.axis=1.25,cex.lab=1.25)
# filled.contour(sigs,rhos,log10(std),xlab=expression(sigma[a]),
               # ylab=expression(sigma[b]),main="log standard deviation",
               # color.palette = topo.colors )
# dev.off()

# pdf("single-variableC.pdf",width=6,height=5)
# par(cex.axis=1.25,cex.lab=1.25)
# filled.contour(sigs,rhos,c,xlab=expression(sigma[a]),
               # ylab=expression(sigma[b]),main="closedness",
               # color.palette = topo.colors )
# dev.off()