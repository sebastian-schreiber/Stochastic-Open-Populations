# The Example-Coral-Matrix.R file for the Methods in Ecology and Evolution (2018) paper
# "The structured demography of open populations in fluctuating environments" 
# by Sebastian Schreiber and Jacob Moore. 
# See the READ-ME.txt for an overview of the files and their dependencies. 

# This file applies the major results of the MEE paper to a matrix population of corals
# using data from Pascual, Mercedes, and Hal Caswell. 
# "The dynamics of a size-classified benthic population with reproductive subsidy." 
# Theoretical Population Biology 39.2 (1991): 129-147.  
# It makes use of the Base-Code.R file

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Define functions
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# load functions to simulate stochastic model (runStochasticOpen), and calculate expected means and covariances (calculateIID)
source("Base-Code.R")
#library(viridis)

# Create function that implements sensitivity calculations (Eqns 9, 10, Appendix S5)
calculateSens_coral = function(out, A, E.A, E.n, dA.high, dA.low, dB.high, dB.low, b, p=freq.good, cov.N, s.high, s.low) {
	
	P.state = c(p, 1-p)				# probability of environmental states
	num.stages = dim(A[,,1])[1]		# number of coral stages
	num.states = length(P.state)		# number of environmental states
	 
	Id = diag(dim(A)[1])
	temp.kron=matrix(0,nrow=num.stages^2,ncol=num.stages^2)
	for(i in 1:num.states){
		temp.kron=temp.kron+P.state[i]*(A[,,i]%x%A[,,i])
	}
	Id.kron = diag(dim(temp.kron)[1])
	 
	# Sensitivity of n* w/respect to s.high and s.low (Eqn 9, Appendix S5)
	dN.high = solve(Id-E.A) %*% (p*(dA.high[,,1]%*%E.n + dB.high[,1]))
    dN.low = solve(Id-E.A) %*% ((1-p)*(dA.low[,,2]%*%E.n + dB.low[,2]))
    
    # Elasticity of n* w/respect to s.high and s.low
    dN.high.elas = dN.high * s.high / E.n
    dN.low.elas = dN.low * s.low / E.n
	 
	 # Sensitivity of cov(n*) w/respect to s.high and s.low (Eqn 10, Appendix S5)
	 dE_AxA.dhigh = p*(dA.high[,,1] %x% A[,,1] + A[,,1] %x% dA.high[,,1])
	 	 	 
	 temp1.high = array(NA, dim=c(9,1,2)) # set up here for first part of eqn 26
	 temp2.high = array(NA, dim=c(9,1,2)) # set up here for second part of eqn 26
	 for(i in 1:num.states){
	 	temp1.high[,,i] = P.state[i]* (  (dA.high[,,i]%*%E.n + dB.high[,i] + (A[,,i]-Id)%*%dN.high) %x% (b[,i]+(A[,,i]-Id)%*%E.n)  )
	 	temp2.high[,,i] = P.state[i]* (  (b[,i]+(A[,,i]-Id)%*%E.n) %x% (dA.high[,,i]%*%E.n + dB.high[,i] + (A[,,i]-Id)%*%dN.high)  )
	 }
	 dCov_Anb.dhigh = temp1.high[,,1]+temp1.high[,,2]+temp2.high[,,1]+temp2.high[,,1]
	 vec.dcovN.high = solve(Id.kron-temp.kron) %*% ( dE_AxA.dhigh %*% as.vector(cov.N) + as.vector(dCov_Anb.dhigh)  )
	 dcov.N.high = matrix(vec.dcovN.high, nrow=num.stages, byrow=T)
	 
	 dE_AxA.dlow = (1-p)*(dA.low[,,2] %x% A[,,2] + A[,,2] %x% dA.low[,,2])	 	 
	 temp1.low = array(NA, dim=c(9,1,2)) # set up here for first part of eqn 26
	 temp2.low = array(NA, dim=c(9,1,2)) # set up here for second part of eqn 26
	 for(i in 1:num.states){
	 	temp1.low[,,i] = P.state[i]* (  (dA.low[,,i]%*%E.n + dB.low[,i] + (A[,,i]-Id)%*%dN.low) %x% (b[,i]+(A[,,i]-Id)%*%E.n)  )
	 	temp2.low[,,i] = P.state[i]* (  (b[,i]+(A[,,i]-Id)%*%E.n) %x% (dA.low[,,i]%*%E.n + dB.low[,i] + (A[,,i]-Id)%*%dN.low)  )
	 }
	 dCov_Anb.dlow = temp1.low[,,1]+temp1.low[,,2]+temp2.low[,,1]+temp2.low[,,1]
	 vec.dcovN.low = solve(Id.kron-temp.kron) %*% ( dE_AxA.dlow %*% as.vector(cov.N) + as.vector(dCov_Anb.dlow)  )
	 dcov.N.low = matrix(vec.dcovN.low, nrow=num.stages, byrow=T)

	 return(list(dN.high=dN.high, dN.low=dN.low, dcov.N.high=dcov.N.high, dcov.N.low=dcov.N.low, dN.high.elas=dN.high.elas, dN.low.elas=dN.low.elas))

}

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set up population model
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# Settlement rates
s.high = 0.04
s.low = 0.01	

# Age-specific area covered
a0 = 5 # cm2
a1 = 30
a2 = 125

# Total area available
A.coral = 120000 # cm2

# Transition matrices
coral = array(NA, dim=c(3,3,2))
	# transitions when good year
	coral[,,1] = matrix(c(-s.high*a0+0.5135, -s.high*a1, -s.high*a2,
				0.2072, 0.6667, 0,
				0,   0.1746, 0.7818), nrow=3, byrow=T)
	# transitions when bad year
	coral[,,2] = matrix(c(-s.low*a0+0.5135, -s.low*a1, -s.low*a2,
				0.2072, 0.6667, 0,
				0,   0.1746, 0.7818), nrow=3, byrow=T)

# recruitment matrices			
coral.input = array(0, dim=c(3,2))
	# recruitment in good year
	coral.input[1,1]=s.high*A.coral
	# recruitment in bad year
	coral.input[1,2]=s.low*A.coral

# define derivative matrices for sensitivity calculations
	# for dA1/ds.high, dA2/ds.low [0 for other cases]
	dA.base = matrix(c(-5, -30, -125,
					0, 0, 0,
					0, 0, 0), nrow=3, byrow=T)
					
	dA.high = array(0, dim=c(3,3,2))
	dA.high[,,1] = dA.base
	dA.low = array(0, dim=c(3,3,2))
	dA.high[,,2] = dA.base

	# for db1/ds.high, db2/ds.low [0 for other cases]
	dB.base = matrix(c(120000,0,0), nrow=3)	
	dB.high = matrix(0, nrow=3, ncol=2)
	dB.high[,1] = dB.base
	dB.low = matrix(0, nrow=3, ncol=2)
	dB.low[,2] = dB.base			
				

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Run model and calculations
# -------------------------------------------------------------------
# -------------------------------------------------------------------

n.initial = matrix(c(1200, 800, 600), ncol=1)	# initial population size
num.steps = 100000								# number of iterations to simulate
freq.good.seq = seq(from=0.05, to=0.95, length=19)	# run model with varying frequencies of a good year

# set up arrays to store data for each iteration
states.all = array(NA, dim=c(num.steps+1, length(freq.good.seq)))
ws.open.all = array(NA, dim=c(3, num.steps+1, length(freq.good.seq)))
n.all = array(NA, dim=c(num.steps+1, length(freq.good.seq)))
lambdaS.all = array(NA, dim=c(length(freq.good.seq)))
E.A.all = array(NA, dim=c(3,3,length(freq.good.seq)))
E.B.all = array(NA, dim=c(3,length(freq.good.seq)))
E.n.all = array(NA, dim=c(3,1,length(freq.good.seq)))
cov.N.all = array(NA, dim=c(3,3, length(freq.good.seq)))
cov.N.all.num = cov.N.all
cor.N.all = cov.N.all
cor.N.all.num = cov.N.all
cov.tau.all = array(NA, dim=c(3,3, length(freq.good.seq)))
cor.tau.all = cov.tau.all
cov.tau.all.num = cov.tau.all
cor.tau.all.num = cov.tau.all
dN.high.all = array(NA, dim=c(3,length(freq.good.seq)))
dN.high.all.elas = dN.high.all
dN.low.all = dN.high.all
dN.low.all.elas = dN.low.all
dcov.N.high.all = array(NA, dim=c(3,3,length(freq.good.seq)))
dcov.N.low.all = dcov.N.high.all

# run through model at different values of p
for(f in 1:length(freq.good.seq)) {
	freq.good = freq.good.seq[f]
	
	# environmental probabilities
	P.coral = c(freq.good, 1-freq.good)
	P.trans.coral = matrix(rep(P.coral, 2), nrow=2, byrow=TRUE)	# because iid
	
	# run simulation
	initial.state = which(P.coral==sample(P.coral, 1, prob=P.coral))[1]
	coral.out = runStochasticOpen(A=coral, B=coral.input, P.state=P.coral, P.trans=P.trans.coral, num.steps, initial.state, n.initial)
	
	tau = 1
	# calculate analytic means, covariances, and correlations
	coral.outIID = calculateIID( A=coral, B=coral.input, P.state=P.coral, tau=tau)

	# calculate covariance and correlation of simulation
	cov.N.all.num[,,f] = cov(t(coral.out$ws.open))
	cor.N.all.num[,,f] = cov2cor(cov.N.all.num[,,f])
	
	# calculate temporal covariance and correlation of simulation
	X=t(coral.out$ws.open)
	X1=X[-((num.steps-tau+2):(num.steps+1)),]
	X2=X[-(1:tau),]
	cov.tau.all.num[,,f] = cov(X1, X2)
	for(i in 1:3) {
		for(j in 1:3) {
			cor.tau.all.num[i,j,f] = cov.tau.all.num[i,j,f] / sqrt(cov.N.all.num[i,i,f]*cov.N.all.num[j,j,f])
		}
	}

	# place output into larger arrays
	states.all[,f] = coral.out$states
	ws.open.all[,,f] = coral.out$ws.open
	n.all[,f] = coral.out$n.open
	lambdaS.all[f] = coral.out$lambdaS
	E.A.all[,,f] = coral.outIID$E.A
	E.B.all[,f] = coral.outIID$E.B
	E.n.all[,,f] = coral.outIID$E.n
	cov.N.all[,,f] = coral.outIID$cov.N
	cor.N.all[,,f] = cov2cor(coral.outIID$cov.N)
	cov.tau.all[,,f] = coral.outIID$cov.tau
	for(i in 1:3) {
		for(j in 1:3) {
			cor.tau.all[i,j,f] = cov.tau.all[i,j,f] / sqrt(cov.N.all[i,i,f]*cov.N.all[j,j,f])
		}
	}

	# calculate sensitivities 
	coral.sens = calculateSens_coral(out=coral.outIID, A= coral, E.A=coral.outIID$E.A, E.n=coral.outIID$E.n, dA.high, dA.low, dB.high, dB.low, b=coral.input, p=freq.good, cov.N=coral.outIID$cov.N, s.high, s.low)

	# store sensitivity calculations
	dN.high.all[,f] = coral.sens$dN.high
	dN.low.all[,f] = coral.sens$dN.low
	dcov.N.high.all[,,f] = coral.sens$dcov.N.high	
	dcov.N.low.all[,,f] = coral.sens$dcov.N.low
	dN.high.all.elas[,f] = coral.sens$dN.high.elas
	dN.low.all.elas[,f] = coral.sens$dN.low.elas

}

# last row of ws.open when running at freq.good of 0 or 1
state.dist.freq.good0 = c(892.3238, 554.7240, 443.8809)
state.dist.freq.good1 = c(1224.5794, 761.2747, 609.1593)

# # -------------------------------------------------------------------
# # -------------------------------------------------------------------
# # Plot output
# # -------------------------------------------------------------------
# # -------------------------------------------------------------------

# misc formatting
ltys=c(1,1,1,2,2,2)
#cols=c(1:3,1:3)
size1=expression(paste('Size class: 0-10 ',cm^2))
size2=expression(paste('Size class: 10-50 ',cm^2))
size3=expression(paste('Size class: 50-200 ',cm^2))

# # for B+W
col.scheme.3 = c('#000000','#737373','#bdbdbd')
col.hist = c('black', 'black', 'black')
dens.val = 0

pdf("coral-dist.pdf",width=7,height=5)
#quartz(width=7,height=5)
par(mfrow=c(2,3),cex.axis=1.25,cex.lab=1.25,mar=c(5.25,4.5,1,1))
i1=1
i2=10
a1=min(ws.open.all[1,,i1],ws.open.all[1,,i2])
b1=max(ws.open.all[1,,i1],ws.open.all[1,,i2])
a2=min(ws.open.all[2,,i1],ws.open.all[2,,i2])
b2=max(ws.open.all[2,,i1],ws.open.all[2,,i2])
a3=min(ws.open.all[3,,i1],ws.open.all[3,,i2])
b3=max(ws.open.all[3,,i1],ws.open.all[3,,i2])
hist(ws.open.all[1,,i1], ylab='frequency', xlab='', main=size1,xlim=c(a1,b1), col=col.hist[1], density=dens.val, freq=FALSE)
abline(v=state.dist.freq.good0[1], col='black', lty=2)
abline(v=state.dist.freq.good1[1], col='black', lty=3)
legend("center","5% good\n years",bty="n", adj=c(-0.4, 0))
hist(ws.open.all[2,,i1], ylab='', xlab='', main=size2,xlim=c(a2,b2), col= col.hist[2], density=dens.val, freq=FALSE)
abline(v=state.dist.freq.good0[2], col='black', lty=2)
abline(v=state.dist.freq.good1[2], col='black', lty=3)
hist(ws.open.all[3,,i1], ylab='', xlab='', main=size3,xlim=c(a3,b3), col= col.hist[3], density=dens.val, freq=FALSE)
abline(v=state.dist.freq.good0[3], col='black', lty=2)
abline(v=state.dist.freq.good1[3], col='black', lty=3)
hist(ws.open.all[1,,i2], ylab='frequency', xlab='population density',main="",xlim=c(a1,b1), col= col.hist[1], density=dens.val, freq=FALSE)
abline(v=state.dist.freq.good0[1], col='black', lty=2)
abline(v=state.dist.freq.good1[1], col='black', lty=3)
legend("right","50% good\n years",bty="n")
hist(ws.open.all[2,,i2], ylab='', xlab='population density',main="",xlim=c(a2,b2), col= col.hist[2], density=dens.val, freq=FALSE)
abline(v=state.dist.freq.good0[2], col='black', lty=2)
abline(v=state.dist.freq.good1[2], col='black', lty=3)
hist(ws.open.all[3,,i2], ylab='', xlab='population density',main="",xlim=c(a3,b3), col= col.hist[3], density=dens.val, freq=FALSE)
abline(v=state.dist.freq.good0[3], col='black', lty=2)
abline(v=state.dist.freq.good1[3], col='black', lty=3)
dev.off()


# Figure of means, sds, correlations, with comparisons of analytic and numeric (Fig 3)
pdf("coral-stats.pdf",width=8,height=6)
#quartz(width=8, height=6)
par(mfrow=c(2,2),cex.axis=1.25,cex.lab=1.25,mar=c(4.5,4.5,1,1))

# equilibria	(Fig 3a)
ltys=c(1,1,1,2,2,2)
y.sims=cbind(colMeans(ws.open.all[1,,]),colMeans(ws.open.all[2,,]),colMeans(ws.open.all[3,,]))
y.analytic=t(E.n.all[1:3,1,])
y=cbind(y.sims,y.analytic)
matplot(freq.good.seq,y.analytic,type="l",lty=ltys,col=col.scheme.3,lwd=2, xlab='frequency of good years', ylab='mean density')
temp=seq(1,length(freq.good.seq),by=3)
matplot(freq.good.seq[temp],y.sims[temp,],type="p",col=col.scheme.3,add=TRUE,pch=4)
legend('top',legend=c(size1,size2,size3), lty=ltys,lwd=2,col=col.scheme.3, cex=0.95,bty="n", inset=c(0, 0.3))
legend('topleft','A',bty="n",cex=1.25)

# variances		(Fig 3b)
plot.start=1
plot.end = length(freq.good.seq)
y=cbind(cov.N.all.num[1,1,plot.start:plot.end], cov.N.all.num[2,2,plot.start:plot.end], cov.N.all.num[3,3,plot.start:plot.end], cov.N.all[1,1,plot.start:plot.end], cov.N.all[2,2,plot.start:plot.end], cov.N.all[3,3,plot.start:plot.end])
ltys=c(1,1,1,2,2,2)
#cols=c(1:3,1:3)
matplot(freq.good.seq,sqrt(y[,4:6]),type="l",lty=ltys,col=col.scheme.3,lwd=2, xlab='frequency of good years', ylab='standard deviation')
temp=seq(1,length(freq.good.seq),by=3)
matplot(freq.good.seq[temp],sqrt(y[temp,1:3]),type="p",col=col.scheme.3,add=TRUE,pch=4)
legend('topleft','B',bty="n",cex=1.25)

# correlations	(Fig 3c)
y=cbind(cor.N.all.num[1,2,plot.start:plot.end], cor.N.all.num[2,3,plot.start:plot.end], cor.N.all.num[1,3,plot.start:plot.end], 
        cor.N.all[1,2,plot.start:plot.end], cor.N.all[2,3,plot.start:plot.end], 
        cor.N.all[1,3,plot.start:plot.end])
ltys=c(6,6,6,6,6,6)
#ltys=c(1,1,1,2,2,2)
#cols=c(1:3,1:3)
matplot(freq.good.seq,y[,4:6],type="l",lty=ltys,col=col.scheme.3,lwd=2, xlab='frequency of good years', ylab='correlation')
temp=seq(1,length(freq.good.seq),by=3)
matplot(freq.good.seq[temp],y[temp,1:3],type="p",col=col.scheme.3,add=TRUE,pch=4)
legend('topleft','C',bty="n",cex=1.25)
size12=expression(paste('0-10 ',' & ','10-50 ',cm^2))
size23=expression(paste('10-50 ',' & ','50-200 ',cm^2))
size31=expression(paste('50-200 ',' & ','0-10 ',cm^2))
legend('center',legend=c(size12,size23,size31), lty=ltys,lwd=2,col=col.scheme.3, cex=0.95,bty="n")

# temporal autocorrelations		(Fig 3d)
y=cbind(cor.tau.all.num[1,1,plot.start:plot.end], cor.tau.all.num[2,2,plot.start:plot.end], cor.tau.all.num[3,3,plot.start:plot.end], cor.tau.all[1,1,plot.start:plot.end], cor.tau.all[2,2,plot.start:plot.end], cor.tau.all[3,3,plot.start:plot.end])
ltys=c(1,1,1,2,2,2)
#cols=c(1:3,1:3)
matplot(freq.good.seq,y[,4:6],type="l",lty=ltys,col=col.scheme.3,lwd=2, 
        xlab='frequency of good years',
        ylab=expression(paste('autocorrelation with ',tau==1)))
legend('topleft','D',bty="n",cex=1.25)
temp=seq(1,length(freq.good.seq),by=3)
matplot(freq.good.seq[temp],y[temp,1:3],type="p",col=col.scheme.3,add=TRUE,pch=4)
dev.off()


# Figure of elasticitiy of mean density and elasticity of standard deviation of mean density (Fig 4)
pdf("coral-elasticity.pdf",width=7,height=3)
par(mfrow=c(1,2),cex.axis=1.25,cex.lab=1.25,mar=c(4.5,4.5,1,1))
ltys=c(1,1,1,2,2,2)
cols=c(1:3,1:3)
size1=expression(paste('0-10 ',cm^2))
size2=expression(paste('10-50 ',cm^2))
size3=expression(paste('50-200 ',cm^2))

# equilibrium	(Fig 4a)
y=cbind(dN.high.all.elas[1,], dN.high.all.elas[2,], dN.high.all.elas[3,], dN.low.all.elas[1,], dN.low.all.elas[2,], dN.low.all.elas[3,])
matplot(freq.good.seq,y,type="l",lty=ltys,col=c('black'),lwd=2, xlab='frequency of good years', ylab='elasticity of mean density')
legend('topright', legend=c("high recruitment","low recruitment"), lty=ltys[c(1,4)], col=1, cex=1,bty="n")


# standard deviation	(Fig 4b)
y=cbind(dcov.N.high.all[1,1,]*s.high, dcov.N.high.all[2,2,]*s.high, dcov.N.high.all[3,3,]*s.high, dcov.N.low.all[1,1,]*s.low, dcov.N.low.all[2,2,]*s.low, dcov.N.low.all[3,3,]*s.low)
y2=cbind(cov.N.all[1,1,plot.start:plot.end], cov.N.all[2,2,plot.start:plot.end], cov.N.all[3,3,plot.start:plot.end])
y2=cbind(y2,y2)
y=y/(2*y2) 
matplot(freq.good.seq,y,type="l",lty=ltys,col=c(col.scheme.3, 'black','black','black'),lwd=2, xlab='frequency of good years', ylab='elasticity of SD')
legend('center', legend=c(size1,size2,size3), lty=ltys, col=col.scheme.3, cex=1,bty="n")

dev.off()

