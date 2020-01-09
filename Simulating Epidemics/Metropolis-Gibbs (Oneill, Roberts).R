#'
#'
#' Metropolis within Gibbs Algorithm 
#' Based on O'neill, Roberts paper '98
#' Uniform proposal for infection times 
#'


SIR.MetropolisGibbs2 = function(N, R, theta.gamma, theta.beta, theta, no.its){
  
  # == Initialise ==
  
  # Infection and Removal rates 
  #beta = rexp(1, rate = theta.beta[1], shape = theta.beta[2])
  #gamma = rgamma(1, rate = theta.gamma[1], shape = theta.gamma[2])
  
  beta = rexp(1, rate = theta.beta)
  gamma = rexp(1, rate = theta.gamma)
  
  # First infection time
  
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
  Q[infected] = rexp(n_I, rate = gamma)
  
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
      Q = Q*1.1
      print(mean(Q))
      # Calculate new infection times
      I = R - Q
    }
  }
  
  # Initial Infected
  k = which.min(I)
  
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
      new = log(beta*sum(I[k] > I & I[k] < R, na.rm = TRUE))  
      logL.list = c(logL.list, new)
      logL = logL + new
    }
  }
  
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
    
    # ==== Update Initial infection time ====
    
    Qk = rexp(1, rate = theta + beta*(N-1) + gamma)
    
    I[k] = min(I[-k]) - Qk
    
    # ==== Updating Infection time(s) ====
    
    # Choose an infected individual at random
    proposed.infected = sample(infected, 1)
    
    # Propose an new infection time for chosen individual
    # by sampling uniformly on (I_k, T)
    I.prop = I
    I.prop[proposed.infected] = runif(1, min = min(I), max = max(R[R < 10000]))
    
    # == Calculate the likelihood for current and proposed infection times ==

    curr.llh = SIR.likelihood(I, R, beta, gamma, log = TRUE)
    prop.llh = SIR.likelihood(I.prop, R, beta, gamma, log = TRUE)
    
    log.u = log(runif(1))
    log.a = (prop.llh$llh) - (curr.llh$llh)
    
    if(log.u < log.a){
      I = I.prop
      logL = prop.llh$logL
      inf.pressure = prop.llh$inf.pressure
      accept = accept + 1
    }
    else{
      logL = curr.llh$logL
      inf.pressure = curr.llh$inf.pressure
    }
    draws[i+1,] = c(beta, gamma, I)
    print(accept)
  }
  accept.rate = accept/no.its
  
  return(list(draws = draws, accept.rate = accept.rate))
}