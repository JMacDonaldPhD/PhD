#'
#'
#'
#' Investigation: Large Epidemics vs. small Epidemics (Popn 200 vs. Popn 20)
#'

# ==== Preamble ====

# Simulated Epidemic Datasets
source("/home/macdona8/Documents/PhD/Simulating Epidemics/Simulation for MWG investigation.R")

# MCMC Algorithms 
source("/home/macdona8/Documents/PhD/Simulating Epidemics/Metropolis-Gibbs (Neal, Roberts).R")
source("/home/macdona8/Documents/PhD/Simulating Epidemics/NCP Epidemic Sampler (Neal, Roberts).R")
source("/home/macdona8/Documents/PhD/Simulating Epidemics/ASIS Sampler.R")

# ==== Centered Parameterisation ====

#' Look at efficiency of Centered Algorithm on Low, Medium and High scale epidemics
#' Use a range of infection time proposals
#' Run until convergence then look at Effective Sample Size per second as a measurement 
#' of efficiency.
#' Start algorithm at true values of beta and gamma used in simulation 

# == Small Scale ==

N = MedPopn.Epidemic$N
R = MedPopn.Epidemic$removal.times

gamma0 = 0.15
beta0 = 0.001
gamma.prior = 'exp' 
theta.gamma = 0.001
beta.prior = 'exp' 
theta.beta = 0.001
rinf.dist = 'exp'
no.its = 10000
burn.in = 1000
#no.proposals = round(c(1/LowPopn.Epidemic$final.size, 0.25, 0.75, 1)*LowPopn.Epidemic$final.size)
no.proposals = rep(1:LowPopn.Epidemic$final.size, 3)
ESS.sec.beta = c()
ESS.sec.gamma = c()
mean.beta = c()
mean.gamma = c()
accept.rate = c()
R0 = c()

run = SIR.MetropolisGibbs(N, R, gamma0, beta0, gamma.prior, theta.gamma, beta.prior, theta.beta, rinf.dist, alpha = 1, 
                    no.its = 20000, burn.in = 1000, no.proposals = 5)
plot(run$draws[,1], type = 'l')
mean(run$draws[,1])
mean(run$draws[,2])
run2 = NCP.Sampler(N, R, gamma0, beta0, gamma.prior, theta.gamma, beta.prior, theta.beta, rinf.dist,
                   no.its = 10000, burn.in = 1000, no.proposals = 5, lambda = 0.3)
mean(run2$draws[,1])
mean(run2$draws[,2])
20*mean(run$draws[,1])/mean(run$draws[,2])


for(i in 1:length(no.proposals)){
  run = SIR.MetropolisGibbs(N, R, gamma0, beta0, gamma.prior, theta.gamma, beta.prior, theta.beta, rinf.dist, alpha = 1, 
                            no.its, burn.in, no.proposals = no.proposals[i])
  ESS.sec.beta[i] = run$ESS.sec[1]
  ESS.sec.gamma[i] = run$ESS.sec[2]
  mean.beta[i] = mean(run$draws[,1])
  mean.gamma[i] = mean(run$draws[,2])
  R0[i] = N*mean.beta/mean.gamma
  accept.rate[i] = run$accept.rate
  print(i)
}


results = data.frame(no.proposals, mean.beta, mean.gamma, R0, accept.rate)
R0 = N*mean.beta/mean.gamma
par(mfrow = c(1,2))

plot(no.proposals, R0, ylim = c(1,2))
abline(a = N*beta0/gamma0, b = 0)
plot(no.proposals, ESS.sec.beta)
plot(no.proposals, ESS.sec.gamma)

# == Medium Scale ==

N = MedPopn.Epidemic$N
R = MedPopn.Epidemic$removal.times
gamma0 = 0.15
beta0 = 0.01
gamma.prior = 'exp' 
theta.gamma = 0.001
beta.prior = 'exp' 
theta.beta = 0.001
rinf.dist = 'exp'
no.its = 10000
burn.in = 1000
#no.proposals = round(c(1/LowPopn.Epidemic$final.size, 0.25, 0.75, 1)*LowPopn.Epidemic$final.size)
no.proposals = round(seq(1:LowPopn.Epidemic, length = 10))
no.proposals = rep(no.proposals, 3)
lambda.values =
ESS.sec.beta = c()
ESS.sec.gamma = c()
mean.beta = c()
mean.gamma = c()
accept.rate = c()

for(i in 1:length(no.proposals)){
  accept.rate = 0
  while(accept.rate < 20 & accept.rate > 40){
    run = NCP(N, R, gamma0 = 0.15, beta0 = 0.001, gamma.prior, theta.gamma, beta.prior, theta.beta, rinf.dist, alpha = 1, 
                              no.its, burn.in, no.proposals = no.proposals[i])
    accept.rate = run$accept.rate
    if(accept.rate )
  }

  ESS.sec.beta[i] = run$ESS.sec[1]
  ESS.sec.gamma[i] = run$ESS.sec[2]
  mean.beta[i] = mean(run$draws[,1])
  mean.gamma[i] = mean(run$draws[,2])
  accept.rate[i] = run$accept.rate
}

results = data.frame(no.proposals, mean.beta, mean.gamma, R0, accept.rate)

par(mfrow = c(1,2))
plot(no.proposals, ESS.sec.beta)
plot(no.proposals, ESS.sec.gamma)


# == Large Scale ==



# ==== Non-Centered Parameterisation ====

# == Small Scale ==
N = LowPopn.Epidemic$N
R = LowPopn.Epidemic$removal.times
gamma0 = 0.15
beta0 = 0.01
gamma.prior = 'exp' 
theta.gamma = 0.001
beta.prior = 'exp' 
theta.beta = 0.001
rinf.dist = 'exp'
no.its = 10000
burn.in = 1000
#no.proposals = round(c(1/LowPopn.Epidemic$final.size, 0.25, 0.75, 1)*LowPopn.Epidemic$final.size)
no.propsals = round(seq(1:LowPopn.Epidemic$final.size, length = 10))
no.proposals = rep(1:LowPopn.Epidemic$final.size, 3)
ESS.sec.beta = c()
ESS.sec.gamma = c()
mean.beta = c()
mean.gamma = c()
accept.rate = c()
R0 = c()
for(i in 1:length(no.proposals)){
  run = NCP.Sampler(N, R, gamma0 = 0.15, beta0 = 0.01, gamma.prior, theta.gamma, beta.prior, theta.beta, rinf.dist, alpha = 1, 
                            no.its, burn.in, no.proposals = no.proposals[i], lambda)
  ESS.sec.beta[i] = run$ESS.sec[1]
  ESS.sec.gamma[i] = run$ESS.sec[2]
  mean.beta[i] = mean(run$draws[,1])
  mean.gamma[i] = mean(run$draws[,2])
  R0[i] = N*mean.beta/mean.gamma
  accept.rate[i] = run$accept.rate
  print(i)
}

results = data.frame(no.proposals, mean.beta, mean.gamma, R0, accept.rate)

par(mfrow = c(1,2))
plot(no.proposals, ESS.sec.beta)
plot(no.proposals, ESS.sec.gamma)

# == Medium Scale ==

# == Large Scale ==



# ==== Partial Non-centering ====










