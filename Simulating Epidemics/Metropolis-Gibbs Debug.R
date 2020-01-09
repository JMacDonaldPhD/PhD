#'
#'
#' Gibbs sampler using Uniform proposal.
#'
#'
#'
#'
#'

# ==== Preamble ====

# Epidemic Simulator
source("/home/macdona8/Documents/PhD/Simulating Epidemics/GSE Simulations.R")

# Functions for Metropolis-Within-Gibbs 
source("/home/macdona8/Documents/PhD/Simulating Epidemics/Functions for MWG.R")


# ==== Simulating Test Data ====
# Set Seed to ensure we get an epidemic which takes off
set.seed(1011)


Sim = GSE.simulation(N = 200, a = 1, gamma = 3, beta = 0.05)

# Extract output from simulation for use in MWG sampler
N = Sim$N
R = Sim$removal.times
I = Sim$infection.times

# Distributions of the various parameters/random variables
beta.prior = "exp" # Infectious Process Prior
gamma.prior = "exp" # Removal rate prior
rinf.dist = "exp" # Length of Infectious Period

# Prior hyperparameters for rate parameters
theta.gamma = 0.001 # rate for prior gamma ~ exp()
theta.beta = 0.001 # rate for prior beta ~ exp()

# Extra Parameters for Infectious Period Distribution
# alpha = 2 # Shape parameter for Gamma/Weibull Distribution

no.its = 3000
# set.seed(35)



# How many Infection times to propose?
no.proposals = 1

# ==== Metropolis-Within-Gibbs Sampler ====

set.seed(sample(1:1000, size = 1))

Start = Sys.time()

# == Initialise ==

# Infection and Removal rates 
#beta = rexp(1, rate = theta.beta[1], shape = theta.beta[2])
#gamma = rgamma(1, rate = theta.gamma[1], shape = theta.gamma[2])

beta = rexp(1, rate = theta.beta)
gamma = rexp(1, rate = theta.gamma)


#return(list(gamma, beta))


# Number of people who were removed 
n_R = sum(R < 10000) 

# For finished epidemics, the number of people who were infected
# is equal to the number who were removed.
#n_I = sum(I < 10000)
n_I = n_R

# Which individuals got infected?
infected = which(R < 10000)

# Rescale so that first removal time is at t = 0
# This time is when the epidemic would 
# begin to be observed, in reality
R[infected] = R[infected] - min(R[infected])

Q = rep(0, N)
I = rep(10000, N)

# Draw infectious periods
Q[infected] = rinf.period(n = n_I, dist = rinf.dist, 
                          theta = switch(rinf.dist, 
                                         exp = gamma, gamma = c(gamma, alpha),
                                         weibull = c(gamma, alpha)))

# Use these infectious periods and the known removal times to 
# to calculate infection times.
I = R - Q

# == Checking validity of drawn infection times ==

# Each infection time must lie within the infectious period of another 
# individual (excluding the initial infected individual)
valid.infected = 0
while(valid.infected < n_I - 1){
  valid.infected = 0
  # Valid infection times check
  for(i in infected){
    if(sum(I[-i] < I[i] & I[i] < R[-i]) > 0){
      valid.infected = valid.infected + 1
    } 
  }
  print(valid.infected)
  # If the infection times do not create a valid epidemic, we can
  # widen the length of the infectious periods so that there is
  # a higher chance of crossover in infectious periods
  if(valid.infected < n_I - 1){
    # Widen infectious periods
    Q = Q*1.05
    # Calculate new infection times
    I = R - Q
  }
}

# Calculate components of posterior parameters for beta and the likelihood
# associated with these infection times and removal times
inf.pressure = 0
logL = 0
logL.list = c()
for(k in which(I<10000)){
  x = 0
  for(j in 1:N){
    x = x + (min(R[k], I[j]) - min(I[j],I[k]))
  }
  inf.pressure = inf.pressure + x
  if(k != which.min(I)){
    new = log(beta*sum(I[k] > I & I[k] < R))  
    logL.list = c(logL.list, new)
    logL = logL + new
  }
}
# , na.rm = TRUE
# Empty matrix to store samples
draws = matrix(NA, nrow = no.its + 1, ncol = N + 2)


# First entry is the intial draws
draws[1,] = c(beta, gamma, I)


# Acceptance Counter
accept = 0
# == The Sampler ==
for(i in 1:no.its){
  
  # == Drawing Beta and Gamma (Gibbs Step) ==
  
  # Update beta and gamma by drawing from their conditional posteriors
  
  # Gamma Prior ---> Gamma Posterior
  #gamma =  rgamma(1, shape = theta.gamma[1] + n_R, rate = theta.gamma[2] + sum(R - I, na.rm = TRUE) )
  #beta  =  rgamma(1, shape = theta.beta[1] + n_I - 1, rate = theta.beta[2] + inf.pressure )
  
  # Exponential Prior ---> Gamma Posterior
  gamma =  rgamma(1, shape = n_R, rate = theta.gamma + sum(R[infected] - I[infected]))
  beta  =  rgamma(1, shape = n_I - 1, rate = theta.beta + inf.pressure)
  #print(c(gamma, sum(R-I, na.rm = TRUE)))
  
  
  # ==== Updating Infection time(s) ====
  
  # Choose an infected individual at random
  proposed.infected = sample(infected, no.proposals)
  
  # Propose a new infectious period
  Q.proposal = rinf.period(no.proposals, dist = rinf.dist, 
                           theta = switch(rinf.dist, 
                                          exp = gamma, gamma = c(gamma, alpha),
                                          weibull = c(gamma, alpha)))
  
  # Use this to calculate the new infection time
  I.prop = I
  I.prop[proposed.infected] = R[proposed.infected] - Q.proposal
  
  
  
  # == Recalculate likelihood components with new beta/gamma values ==
  inf.pressure = 0
  logL = 0
  logL.list = c()
  for(k in which(I<10000)){
    x = 0
    for(j in 1:N){
      x = x + (min(R[k], I[j]) - min(I[j],I[k]))
    }
    inf.pressure = inf.pressure + x
    if(k != which.min(I)){
      new = log(beta*sum(I[k] > I & I[k] < R))  
      logL.list = c(logL.list, new)
      logL = logL + new
    }
  }
  #, na.rm = TRUE
  # == Calculate the likelihood components with proposed infection time(s) ==
  inf.pressure.prop = 0 
  logL.prop =  0
  logL.prop.list = c()
  for(k in which(I.prop<10000)){
    x = 0
    for(j in 1:N){
      x = x + (min(R[k], I.prop[j]) - min(I.prop[j],I.prop[k]))
    }
    inf.pressure.prop = inf.pressure.prop + x
    if(k != which.min(I.prop)){
      new = log(beta*sum(I.prop[k] > I.prop & I.prop[k] < R))
      logL.prop.list = c(logL.prop.list, new) 
      logL.prop =  logL.prop + new
    }
  }
  # , na.rm = TRUE
  log.u = log(runif(1))
  
  # Log M-H acceptance probability 
  log.a = (logL.prop - beta*inf.pressure.prop) - (logL - beta*inf.pressure)
  
  if(log.u < log.a){
    I = I.prop
    logL = logL.prop
    inf.pressure = inf.pressure.prop
    accept = accept + 1
  }
  
  draws[i+1,] = c(beta, gamma, I)
  print(accept)

}


End = Sys.time()
time.taken = as.numeric(End - Start)
ESS = effectiveSize(draws[,1:2])
ESS.sec = ESS/time.taken 
accept.rate = accept/no.its



# ==== Messing Area ====

#plot(draws[-(1:100),2], type ='l')
#acf(draws[1:10,2])