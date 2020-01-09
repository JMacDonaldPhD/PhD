#'
#' PhD
#' SIR Epidemic Likelihood
#'
#'

# beta := infection rate
# gamma :=  removal rate
# I := infection times
# R := removal times

SIR.likelihood = function(I, R, beta, gamma, log = TRUE){
  N = length(I)
  n_I = length(I) - sum(I == 10000)# Number of infected 
  n_R = length(R) - sum(R == 10000)# Number of removed 
  # Calculate the infectious pressure felt by each individual who is infected, at 
  # the time of their infection. Not including first individual
  #for(i in which(I<10000)[-min(I)]){
  #  infection.pressure = c(infection.pressure, beta*sum(I[i] > I & I[i] < R))
  #}
  # Calculate the first part of the likelihood, consiting of the contribution
  # of infections and removals.
  #L1 = prod(infection.pressure)*gamma^(n_R)
  
  #for(i in which(I<10000)){
  #  x = 0
  #  for(j in 1:N){
  #    x = x + (min(R[i], I[j]) - min(I[j],I[i]))
  #  }
  #  infectious.pressure = c(infectious.pressure, x)
  #}
  
  inf.pressure = 0
  logL = 0
  for(k in which(I<10000)){
    x = 0
    for(j in 1:N){
      x = x + (min(R[k], I[j]) - min(I[j],I[k]))
    }
    inf.pressure = inf.pressure + x
    if(k != which.min(I)){
      new = log(beta*sum(I[k] > I & I[k] < R, na.rm = TRUE))  
      logL = logL + new
    }
  }
  
  # == Combine components of the likelihood == 
  gamma.llh = -gamma*sum(R-I) + n_R*log(gamma)
  beta.llh = logL - beta*inf.pressure
  llh = gamma.llh + beta.llh
  
  
  return(list(llh = llh, logL = logL, inf.pressure = inf.pressure))  
}

