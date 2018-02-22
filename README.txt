The files in this repository are the R code used to conduct the analysis and generate the figures for the Methods in Ecology and Evolution (2018) paper "The structured demography of open populations in fluctuating environments" by Sebastian Schreiber and Jacob Moore. Below is a brief description of the files and their interdependency. 

Base-Code.R This file defines two main functions, runStochasticOpen and calculateIID which are used, in most of the other R files. 

	run StochasticOpen  simulates stochastic, structured population models of the types described in the paper and returns the environmental states, population states, and an estimate of the Lyapunov exponent. 

	calculateIID calculates the means and covariances of the stationary distribution of the models for uncorrelated, stationary environments.


Example-Coral-Matrix.R  This file applies the major results of the MEE paper to a matrix population of corals using data from Pascual, Mercedes, and Hal Caswell.  "The dynamics of a size-classified benthic population with reproductive subsidy." Theoretical Population Biology 39.2 (1991): 129-147. It makes use of the Base-Code.R file.

Example-Clam-IPM.R  This file applies the major results of the MEE paper to an integral projection model (IPM)  of giant clams, using data from Yau, Annie J., Hunter S. Lenihan, and Bruce E. Kendall. "Fishery management priorities vary with self‐recruitment in sedentary marine populations." Ecological Applications 24.6 (2014): 1490-1504. It makes use of the Base-Code.R file.


Example-SingleVariable.R This file studies the simple, scalar model n(t+1)=exp(mua+sigmaa*Z1(t+1))*n(t)+exp(mub+sigmab*Z2(t+1)) where Z1, Z2 are standard normals. It makes no use of the Base-Code.R file. 

