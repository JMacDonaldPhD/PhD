#'
#' PhD
#' GSE Simulation
#'
#'


# Population of N individuals
# Start with 1 infected individual, with ID one
# Simulate infectious contact time with each other individuals for every individual.
# However, not everyone will have an infectious period, this infectious contact time is
# required if an individual does have an infectious period
# Infectious contact between two individuals should be simulated for the infectious
# period of the infectious person. Therefore, every individual will have an infectious
# contact wirh every other individual twice but only the first one will be important.
# Of course if neither person gets infected, in reality they will not have infectious
# contact at all.
# One infectious contact time relates to the infectious period of
# one individual, the other relates to the infectious period of the other individual.
# W ~ Exp(beta)


#set.seed(100027)
#Sim = GSE.simulation(N = 200, a = 1, gamma = 3, beta = 0.025)
GSE.simulation = function(N, a, gamma, beta){

  # This status vector will keep track of which set
  # each individual is currently in.
  # Susceptible == 0
  # Infectious == 1
  # Removed == 2
  status = c(rep(1, a), rep(0, N - a))

  # Inf.times will store the infection time for each
  # indiviual. If an individual as not yet been infected
  # there infectious time will be denoted by
  inf.times = c(rep(0, a), rep(NA, N - a))

  # Simulate infectious period lengths foreach
  # individual in the population
  Q = rexp(N, rate = gamma)

  # Simulate contact times
  W = matrix(rexp(N^2, rate = beta), nrow = N, ncol = N)
  diag(W) = rep(10000, N)

  W_2 = matrix(NA, nrow  = N, ncol = N)

  for(i in 1:N){
    for(j in 1:N){
      if(W[i,j] < Q[i]){
        W_2[i,j] = W[i,j]
      }
      else{
        W_2[i,j] = 10000
      }
    }
  }
  # Set final size to zero
  # As infections occure this counter will increase, with the
  # final size at the end of the epidemic, representing the
  # number of people who were infected.
  final.size = 0

  if(min(W_2[1,]) == 10000){
    return(list(final.size = final.size, removal.times = c(Q[1:a], rep(NA, N - a)),
           infection.times = c(rep(0, a), rep(NA, N - a))))
  }
  #
  while(sum(status == 0) > 0 & sum(status == 1) > 0 & sum(status == 1) < N){
    # Set of susceptibles
    S = which(status == 0)
    # Set of Infected
    I = which(status == 1)

    # Reduced Infectious Contact matrix
    W_2.S.I = W_2[I, S, drop = F]

    # Infection times matrix
    A_inf = matrix(inf.times[I], nrow = length(I), ncol = length(S), byrow = FALSE)

    # Potential times of next infections
    M = A_inf + W_2.S.I

    next.inf.time = min(M)
    if(next.inf.time < 10000){
      min.index = which(M == min(M), arr.ind = TRUE)
      infector = I[min.index[1]]
      infected = S[min.index[2]]
      print(min.index)
      print(infector)
      print(infected)
      status[infected] = 1
      inf.times[infected] = next.inf.time

      if(min(W_2[infected,]) == 10000){
        status[infected] = 2
      }

      # If the infector has no more infectious contact with susceptible individuals,
      # change status to removed.
      if(ncol(W_2.S.I) > 1){
        if(min(W_2.S.I[min.index[1], -min.index[2]]) == 10000){
          status[infector] = 2
          I = I[-min.index[1]]
        }
      }
      final.size = final.size + 1
    }
    else{
      status[I] = 2
    }
  }
  rem.times = inf.times + Q
  inf.times[-which(inf.times >= 0)] = 10000
  rem.times[-which(rem.times >= 0)] = 10000
  return(list(final.size = final.size, removal.times = rem.times, infection.times = inf.times))
}
