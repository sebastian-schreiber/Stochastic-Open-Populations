# The Base-Code.R file for the Methods in Ecology and Evolution (2018) paper "The structured 
# demography of open populations in fluctuating environments" by Sebastian Schreiber 
# and Jacob Moore. See the READ-ME.txt for an overview of the files and their dependencies. 

# In this file:
#	1) defines a function, runStochasticOpen, to simulate a stochastic, structured population (matrix or IPM). 
#	2) defines a function, calculateIID, to calculate main results from the manuscript for uncorrelated, stationary environments, including Eqn 4 (expected mean), Eqns 5,6 (covariance of expected mean), and Eqn 8 (temporal covariance of expected mean) 

# runStochasticOpen requires:
#	A: transition matrix
# 	B: external inputs
#	P.state: probability of environmental states
# 	P.trans: conditional probability of transitioning between environmental states
#	num.steps: number of time steps to run the simulation
# 	initial.state: the starting environmental state
# 	n.initial: starting vector of population sizes at each stage/size

# runStochasticOpen returns:
#	states: vector of environmental states at each time step
# 	ws.open: matrix of population sizes for each state/size at each time step for an open population
# 	n.open: vector of total population size at each time step for an open population
#	lambdaS: an estimate of the stochastic population growth rate for a closed population


# calculateIID requires:
#	A: transition matrix
# 	B: external inputs
#	P.state: probability of environmental states
# 	tau: time lag for temporal covariance

# calculateIID returns:
#	E.A: expected transition matrix
#	E.B: expected input vector
#	E.n: expected population sizes
#	cov.N: covariance of expected population sizes
# 	cov.tau: temporal covariance of expected population sizes for given tau



# clear workspace
rm (list = ls ())

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Define functions to run simulation and make i.i.d. calculations
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# Create function that simulates the stochastic population
runStochasticOpen = function(A, B, P.state, P.trans, num.steps, initial.state, n.initial) {
	num.states = length(P.state)		# number of environmental states
	num.stages = dim(A[,,1])[1]		# number of population stages / sizes
	
	# define variables
	lyap = numeric(num.steps+1)						# change in population size each time step in closed population
	n.open = numeric(num.steps+1)					# total population each time step for open population
	states = numeric(num.steps+1)					# environmental state at each time step
	ws.open = matrix(0,num.stages,num.steps+1)		# hold population vector at each time step, open population
	ws.closed = matrix(0,num.stages,num.steps+1)	# hold population vector at each time step, closed population
	
	# define initial conditions
	state = initial.state			# initial state
	states[1]=state 				# add to states vector

	w.open = n.initial 		# initial population vector, open population
	ws.open[,1] = w.open	
	n.open[1]=sum(w.open)	# total population size at time=1, open population
	
	w.closed = w.open			# initial population vector, closed population
	lyap[1]=sum(abs(w.closed))
		
	# simulate population 
	for(t in 1:num.steps) {
		# calculate population size for open population, and store
		w.open = A[,,state] %*% w.open + B[,state]
		ws.open[,t+1]=w.open
		n.open[t+1]=sum(w.open)
		
		# calculate change population size in closed population, and store
		temp.closed = A[,,state] %*% w.closed
		w.closed = temp.closed/sum(abs(temp.closed))
		lyap[t+1]=sum(abs(temp.closed))
		
		# update environmental state
		state = sample(c(1:num.states), 1, prob=P.trans[state,]) 
		states[t+1]=state
	}

	# verify that in absence of input, stochastic population growth rate, lambdaS, <0
	lambdaS = mean(log(lyap))

	# return relevant variables
	return(list(states=states, ws.open=ws.open, n.open=n.open, lambdaS=lambdaS))
}

# Create function that calculates main results from paper for uncorrelated environments (Eqns 4, 5, 6, 8)
calculateIID = function(A, B, P.state, tau) {

	num.states = length(P.state)		# number of environmental states
	num.stages = dim(A[,,1])[1]		# number of population stages / sizes

	# Calculate expected transition matrix (E.A), input vector (E.B), and mean (E.n)
	Id = diag(dim(A)[1])
	E.A = 0
	E.B = 0
	for(i in 1:num.states) {
		E.A = E.A + P.state[i]*A[,,i]
		E.B = E.B + P.state[i]*B[,i]
	}
	E.n = solve(Id-E.A) %*% E.B 		# expected population size of each stage; Eqn 4
	
	# Calculate the covariance of n.hat; Eqn 5
	beta.mat = matrix(0, nrow=num.stages, ncol=num.states)	# let beta = An+b
	cov.beta = matrix(0, nrow=num.stages, ncol=num.stages)	# Cov[An+b]
	for(j in 1:num.states) {
		beta.mat[,j] = A[,,j]%*%E.n + B[,j] - E.n
		cov.beta = cov.beta + P.state[j] * (beta.mat[,j]%*%t(beta.mat[,j]))
	}
	temp.kron=matrix(0,nrow=num.stages^2,ncol=num.stages^2)
	for(i in 1:num.states){
		temp.kron=temp.kron+P.state[i]*(A[,,i]%x%A[,,i])
	}
	Id.kron = diag(dim(temp.kron)[1])
	vec.covN = solve(Id.kron-temp.kron) %*% as.vector(cov.beta)		# Eqn 6
	cov.N = matrix(vec.covN, nrow=num.stages, byrow=T)				# Eqn 5
	
	# Calculate the temporal covariance of n.hat; Eqn 8
	cov.tau = cov.N %*% (t(E.A)^tau)
	
	return(list(E.A=E.A, E.B=E.B, E.n=E.n, cov.N=cov.N, cov.tau=cov.tau))	
}


