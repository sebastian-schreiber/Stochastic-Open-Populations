# The Example-Clam-IPM.R file for the Methods in Ecology and Evolution (2018) paper 
# "The structured demography of open populations in fluctuating environments" 
# by Sebastian Schreiber and Jacob Moore. 
# See the READ-ME.txt for an overview of the files and their dependencies. 


# This file applies the major results of the MEE paper to an integral projection model (IPM) 
# of giant clams, using data from 
# Yau, Annie J., Hunter S. Lenihan, and Bruce E. Kendall.
# "Fishery management priorities vary with self‐recruitment in sedentary marine populations." # Ecological Applications 24.6 (2014): 1490-1504.

rm (list = ls ())

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Define functions
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# load functions to simulate stochastic model (runStochasticOpen), and calculate expected means and covariances (calculateIID)
source("Base-Code.R")

# Create function that implements sensitivity calculations (Eqns 9, 10, Appendix S5)
calculateSens_clam = function(A, E.A, b, E.n) {
	
	num.stages = dim(A[,,1])[1]		# number of sizes in the population
	num.states = dim(b)[2]			# number of environmental states
	
	Id = diag(dim(A)[1])
	temp.kron = A[,,1] %x% A[,,1] 		# since all the A's are the same
	Id.kron = diag(dim(temp.kron)[1])
	
	# Sensitivity of n* when increasing pi and decreasing pj!=pi; Eqn 9, Appendix S5
	dN.p1 = solve(Id-E.A) %*% ( b[,1]-(b[,2]+b[,3]+b[,4])/3)
	dN.p2 = solve(Id-E.A) %*% ( b[,2]-(b[,1]+b[,3]+b[,4])/3)
	dN.p3 = solve(Id-E.A) %*% ( b[,3]-(b[,2]+b[,1]+b[,4])/3)
	dN.p4 = solve(Id-E.A) %*% ( b[,4]-(b[,2]+b[,3]+b[,1])/3)
	
	# Sensitivity of log n*
	dN.p1.log = dN.p1 / E.n
	dN.p2.log = dN.p2/ E.n
	dN.p3.log = dN.p3/ E.n
	dN.p4.log = dN.p4 / E.n
	
	# Sensitivity of cov(n*) when increasing pi and decreasing pj!=pi; Eqn 10, Appendix S5
	dcovN.p1.temp= (b[,1] + (E.A-Id)%*%E.n) %x% (b[,1] + (E.A-Id)%*%E.n)	
	dcovN.p2.temp= (b[,2] + (E.A-Id)%*%E.n) %x% (b[,2] + (E.A-Id)%*%E.n)
	dcovN.p3.temp=(b[,3] + (E.A-Id)%*%E.n) %x% (b[,3] + (E.A-Id)%*%E.n)
	dcovN.p4.temp= (b[,4] + (E.A-Id)%*%E.n) %x% (b[,4] + (E.A-Id)%*%E.n)	
	
	dCov_Anb.p1 = dcovN.p1.temp-(dcovN.p2.temp+dcovN.p3.temp+dcovN.p4.temp)/3	
	vec.dcovN.p1 = solve(Id.kron-temp.kron) %*% as.vector(dCov_Anb.p1)
	dcov.N.p1 = matrix(vec.dcovN.p1, nrow=num.stages, byrow=T)
	
	dCov_Anb.p2 = dcovN.p2.temp-(dcovN.p1.temp+dcovN.p3.temp+dcovN.p4.temp)/3		
	vec.dcovN.p2 = solve(Id.kron-temp.kron) %*% as.vector(dCov_Anb.p2)
	dcov.N.p2 = matrix(vec.dcovN.p2, nrow=num.stages, byrow=T)
	
	dCov_Anb.p3 = dcovN.p4.temp-(dcovN.p2.temp+dcovN.p1.temp+dcovN.p4.temp)/3		
	vec.dcovN.p3 = solve(Id.kron-temp.kron) %*% as.vector(dCov_Anb.p3)
	dcov.N.p3 = matrix(vec.dcovN.p3, nrow=num.stages, byrow=T)
	
	dCov_Anb.p4 = dcovN.p4.temp-(dcovN.p3.temp+dcovN.p3.temp+dcovN.p2.temp)/3	
	vec.dcovN.p4 = solve(Id.kron-temp.kron) %*% as.vector(dCov_Anb.p4)
	dcov.N.p4 = matrix(vec.dcovN.p4, nrow=num.stages, byrow=T)
	
	
	return(list(dN.p1=dN.p1, dN.p2=dN.p2, dN.p3=dN.p3, dN.p4=dN.p4, dN.p1.log=dN.p1.log, dN.p2.log=dN.p2.log, dN.p3.log=dN.p4.log, dN.p4.log=dN.p4.log, dcov.N.p1=dcov.N.p1, dcov.N.p2=dcov.N.p2, dcov.N.p3=dcov.N.p3, dcov.N.p4=dcov.N.p4))
	
}

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Define demographic functions for IPM
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# survival function (probability of surviving)
	s.x = function(currentSize) {
		a.s = 1.29
		b.s = -19.14
		u = exp(a.s + (b.s / (currentSize+3.86)))
		return(u/(1+u))  
	}
# growth function
	g.yx = function(futureSize, currentSize) {	
		a.g = 21.15
		b.g = 0.869
		a.v = 68.45
		b.v = -0.275		
		growth.mean = a.g + b.g*currentSize
		growth.var = a.v + b.v*currentSize
		return(dnorm(futureSize, mean=growth.mean, sd=sqrt(growth.var)))
	}	

# fecundity function
	f.yx = function(futureSize, currentSize, cf) {
		a.f = 3.84e-8
		b.f = 3.63
		mu.rec = 28.3
		sigma.rec = 10.6
		
		if(currentSize<66.1) {		# only reproduce if larger than 66.1 mm
			mean.num.offspring=0
		} else {
			mean.num.offspring = a.f*cf*(currentSize^b.f)
		}
		return(mean.num.offspring * dnorm(futureSize, mean=mu.rec, sd=sigma.rec))
	}  

		
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Set up IPM
# -------------------------------------------------------------------
# -------------------------------------------------------------------

min.size = 1 
max.size = 200
deltax = 4		# step size for discretization (main value 2) - for sensitivity also used 32,16,10,8,4
	
xs = seq(min.size, max.size, by=deltax)	# vector of sizes
n = length(xs)-1
xs.mid = 0.5 * (xs[-1] + xs[-(n+1)])		# midpoint of sizes

# conversion factors
LtoB = exp((log(xs.mid/10)+0.9495)/3.1634)
LtoB.mat = LtoB %*% t(LtoB)

# growth and survival kernels
G=deltax*outer(xs.mid,xs.mid,g.yx) # growth kernel
S=s.x(xs.mid) # survival kernel
P=G
for(i in 1:n) P[,i]=G[,i]*S[i]  # growth kernel conditional on survival


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Run stochastic IPM
# -------------------------------------------------------------------
# -------------------------------------------------------------------
	
	# Specify value for cf:
	cf = 0.25
	
    num.steps = 5000 # number of iterations to simulate

	vect.fyx = Vectorize(f.yx)
	F=deltax*outer(xs.mid,xs.mid,vect.fyx, cf) # reproduction kernel
	K.local=P+F # full kernel		
	
	# transition matrices (all the same)
	clam = array(K.local, dim=c(dim(K.local),4))

	# recruitment matrices
	clam.input = array(0, dim=c(dim(K.local)[1], 4))
	recruit.vect = c(84, 167, 178, 362)
	for(i in 1:length(recruit.vect)) {
		clam.input[,i]= recruit.vect[i] * dnorm(xs.mid, mean=28.3, sd=10.6) 
	}

	# environmental probabilities
	P.clam = c(0.25, 0.25, 0.25, 0.25)
    initial.state = which(P.clam==sample(P.clam, 1, prob=P.clam))[1]
	P.trans.clam = matrix(rep(P.clam, 4), nrow=4, byrow=TRUE)	# equal chance of moving to any other environment

	# starting as 'new' population, from deterministic recruit input
	n.initial = matrix(198 * dnorm(xs.mid, mean=28.3, sd=10.6), ncol=1)
	
	# run model
	clam.out = runStochasticOpen(A=clam, B=clam.input, P.state=P.clam, P.trans=P.trans.clam, num.steps, initial.state, n.initial)


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Run calculations at different levels of local retention, cf
# -------------------------------------------------------------------
# -------------------------------------------------------------------

cf.seq = seq(from=0, to=0.78, length=20)	# vary the degree of local retention, cf
	
# set up arrays to store data for each iteration, and each value of cf
num.stages = n 
E.A.all = array(NA, dim=c(num.stages,num.stages,length(cf.seq)))
E.B.all = array(NA, dim=c(num.stages,length(cf.seq)))
E.n.all = array(NA, dim=c(num.stages,1,length(cf.seq)))
cov.N.all = array(NA, dim=c(num.stages,num.stages, length(cf.seq)))
cov.N.all.num = cov.N.all
cor.N.all = cov.N.all
cor.N.all.num = cov.N.all
cov.tau.all = array(NA, dim=c(num.stages,num.stages, length(cf.seq)))
cor.tau.all = cov.tau.all
cov.tau.all.num = cov.tau.all
cor.tau.all.num = cov.tau.all	

dN.p1.all = array(NA, dim=c(num.stages, length(cf.seq)))
dN.p1.log.all = dN.p1.all
dN.p2.all = array(NA, dim=c(num.stages, length(cf.seq)))
dN.p2.log.all = dN.p2.all
dN.p3.all = array(NA, dim=c(num.stages, length(cf.seq)))
dN.p3.log.all = dN.p3.all
dN.p4.all = array(NA, dim=c(num.stages, length(cf.seq)))
dN.p4.log.all = dN.p4.all
dcov.N.p1.all = array(NA, dim=c(num.stages, num.stages, length(cf.seq)))
dcov.N.p2.all = array(NA, dim=c(num.stages, num.stages, length(cf.seq)))
dcov.N.p3.all = array(NA, dim=c(num.stages, num.stages, length(cf.seq)))
dcov.N.p4.all = array(NA, dim=c(num.stages, num.stages, length(cf.seq)))

N.total.all = array(NA, dim=c(length(cf.seq)))
var.N.total.all = array(NA, dim=c(length(cf.seq)))
sd.N.total.all = array(NA, dim=c(length(cf.seq)))
B.total.all = array(NA, dim=c(length(cf.seq)))
var.B.total.all = array(NA, dim=c(length(cf.seq)))
sd.B.total.all = array(NA, dim=c(length(cf.seq)))

# run through model at different values of cf
for(c in 1:length(cf.seq)) {
	cf = cf.seq[c]

	vect.fyx = Vectorize(f.yx)
	F=deltax*outer(xs.mid,xs.mid,vect.fyx, cf) # reproduction kernel
	K.local=P+F # full kernel		
	
	# transition matrices (all the same)
	clam = array(K.local, dim=c(dim(K.local),4))

	# recruitment matrices
	clam.input = array(0, dim=c(dim(K.local)[1], 4))
	recruit.vect = c(84, 167, 178, 362)
	for(i in 1:length(recruit.vect)) {
		clam.input[,i]= recruit.vect[i] * dnorm(xs.mid, mean=28.3, sd=10.6) 
	}

	# environmental probabilities
	P.clam = c(0.25, 0.25, 0.25, 0.25)
	
	tau = 1
	# calculate means and covariances
	clam.outIID = calculateIID( A=clam, B=clam.input, P.state=P.clam, tau=tau)
	
	# store population size
	N.total.all[c] = sum(clam.outIID$E.n)
	var.N.total.all[c] = sum(diag(clam.outIID$cov.N)) 
	sd.N.total.all[c] = sqrt(var.N.total.all[c])
	
	# store biomass
	B.total.all[c] = sum(clam.outIID$E.n * LtoB) 
	var.B.total.all[c] = sum(LtoB.mat * clam.outIID$cov.N)
	sd.B.total.all[c] = sqrt(var.B.total.all[c])
	
	# store calculations
	E.A.all[,,c] = clam.outIID$E.A
	E.B.all[,c] = clam.outIID$E.B
	E.n.all[,,c] = clam.outIID$E.n
	cov.N.all[,,c] = clam.outIID$cov.N
	cor.N.all[,,c] = cov2cor(clam.outIID$cov.N)
	cov.tau.all[,,c] = clam.outIID$cov.tau
	for(i in 1:3) {
		for(j in 1:3) {
			cor.tau.all[i,j,c] = cov.tau.all[i,j,c] / sqrt(cov.N.all[i,i,c]*cov.N.all[j,j,c])
		}
	}
	
	# calculate sensitivities
	clam.sens = calculateSens_clam(A=clam, E.A=clam.outIID$E.A, b=clam.input, E.n=clam.outIID$E.n)
	
	# store sensitivity calculations
		# sensitivity of mean
		dN.p1.all[,c]=clam.sens$dN.p1
		dN.p2.all[,c]=clam.sens$dN.p2
		dN.p3.all[,c]=clam.sens$dN.p3
		dN.p4.all[,c]=clam.sens$dN.p4
		
		# sensitivity of log mean
		dN.p1.log.all[,c]=clam.sens$dN.p2.log
		dN.p2.log.all[,c]=clam.sens$dN.p2.log
		dN.p3.log.all[,c]=clam.sens$dN.p3.log
		dN.p4.log.all[,c]=clam.sens$dN.p4.log
		
		# sensitivity of cov(n*)
		dcov.N.p1.all[,,c]=clam.sens$dcov.N.p1
		dcov.N.p2.all[,,c]=clam.sens$dcov.N.p2
		dcov.N.p3.all[,,c]=clam.sens$dcov.N.p3
		dcov.N.p4.all[,,c]=clam.sens$dcov.N.p4

}

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Save and/or load data
# -------------------------------------------------------------------
# -------------------------------------------------------------------

## to save the output to avoid re-running, use the following save command
#save.image(file="clams-5.Rdata")

# To avoid re-running use the following load command
# load(file="clams.Rdata")

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Plot output
# -------------------------------------------------------------------
# -------------------------------------------------------------------

col.scheme.map = grey.colors(10, start=0.2, end=0.95)
col.scheme.4 = c('#cccccc','#969696','#636363','#252525')
ltys = c(1,1,1,1)

## The following function needs to be included to plot the color scale: 
# # #This function creates a color scale for use with e.g. the image()
# #function. Input parameters should be consistent with those
# #used in the corresponding image plot. The "horiz" argument
# #defines whether the scale is horizonal(=TRUE) or vertical(=FALSE).
# #Depending on the orientation, x- or y-limits may be defined that
# #are different from the z-limits and will reduce the range of
# #colors displayed.
image.scale <- function(z, zlim, col = heat.colors(12),
breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
 if(!missing(breaks)){
  if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
 }
 if(missing(breaks) & !missing(zlim)){
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
 }
 if(missing(breaks) & missing(zlim)){
  zlim <- range(z, na.rm=TRUE)
  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 poly <- vector(mode="list", length(col))
 for(i in seq(poly)){
  poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
 }
 xaxt <- ifelse(horiz, "s", "n")
 yaxt <- ifelse(horiz, "n", "s")
 if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
 if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
 if(missing(xlim)) xlim=XLIM
 if(missing(ylim)) ylim=YLIM
 plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
 for(i in seq(poly)){
  if(horiz){
   polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
  }
  if(!horiz){
   polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
  }
 }
}
 
pdf("clam_stuff.pdf",width=8.5,height=7)
par(cex.lab=1.25,cex.axis=1.25, mar=c(4.1,4.1,3,2))
layout(matrix(c(1,2,3,4,5,6,7,8,9,9,9,9), nrow=3, ncol=4,byrow=TRUE), widths=c(2.5, 2.5, 2.5, 1), heights=c(2,2,3))
#(Fig 5, top panel)
c1 = 1
plot(xs.mid, E.n.all[,1,c1], type='l', lwd=2, ylab='mean density', xlab='size')
legend("topright", expression(paste(c[f], '=0')), bty="n",cex=1.5)
plot(xs.mid,sqrt(diag(cov.N.all[,,c1]))/E.n.all[,1,c1], type='l', lwd=2, ylab='coefficient of variation', xlab='size')	
image(xs.mid, xs.mid, cor.N.all[,,c1], xlab='size (t)', ylab='size (t+1)', col=col.scheme.map, zlim=c(0,1))
image.scale(c(0,1), col=col.scheme.map, horiz=FALSE, ylab='', xlab='', yaxt='n')
axis(2,at=c(0,0.5,1))
#(Fig 5, middle panel)
c3 = 19
plot(xs.mid, E.n.all[,1,c3], type='l', lwd=2, ylab='mean density', xlab='size')
legend("topright", expression(paste(c[f], '=0.7389')), bty="n",cex=1.5)
plot(xs.mid,sqrt(diag(cov.N.all[,,c3]))/E.n.all[,1,c3], type='l', lwd=2, ylab='coefficient of variation', xlab='size')
image(xs.mid, xs.mid, cor.N.all[,,c3], xlab='size (t)', ylab='size (t+1)', col=col.scheme.map, zlim=c(0,1))
image.scale(range(cor.N.all[,,c3]), col=col.scheme.map[6:10], horiz=FALSE, ylab='', xlab='', yaxt='n')
axis(2,at=c(0.562,0.78,1))
#(Fig 5, bottom panel)
plot(cf.seq,(log(B.total.all-sd.B.total.all)), type='l', lty=2, xlab=expression(c[f]), ylab='total log biomass (g)', lwd=1)
lines(cf.seq,log(B.total.all), lty=1, lwd=2)
lines(cf.seq,log(B.total.all+sd.B.total.all), lty=2)
dev.off()
  
  
# Figure of sensitivity of mean density and sensitivity of standard deviation of mean density (Fig 6)
pdf("clam_sensitivities.pdf", width=8, height=4)
#quartz(width=8, height=4)
par(cex.lab=1.25,cex.axis=1.25,mfrow=c(1,2))
c1 = 1
y1 = cbind(dN.p1.all[,c1], dN.p2.all[,c1], dN.p3.all[,c1], dN.p4.all[,c1])
y1.sd = cbind(diag(dcov.N.p1.all[,,c1]), 
              diag(dcov.N.p2.all[,,c1]),
              diag(dcov.N.p3.all[,,c1]), 
              diag(dcov.N.p4.all[,,c1]))

matplot(xs.mid, y1, type='l', lwd=3, ylab='sensitivity of mean density', xlab='size', col=col.scheme.4, lty=ltys)
abline(h=0, lty=2)
matplot(xs.mid, y1.sd, type='l', lwd=3, ylab='sensitivity of SD', xlab='size', col=col.scheme.4, lty=ltys)
abline(h=0, lty=2)
legend('topright', legend=c(expression(paste(p[1], ' (worst year)', sep='')), expression(p[2]), expression(p[3]), expression(paste(p[4], ' (best year)', sep=''))), lwd=c(rep(2, 4)), col=col.scheme.4, lty=c(1,1,1,1),bty="n")
dev.off()
  
  
  
