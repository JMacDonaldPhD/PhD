#'
#' PhD
#' Experiments involving MGs samplers for epidemics
#'
#'

# ==== Preamble ====
# SIR Simulation
source("/home/macdona8/Documents/PhD/Simulating Epidemics/GSE Simulations.R")

# SIR Likelihood function
source("/home/macdona8/Documents/PhD/Simulating Epidemics/SIR Likelihood.R")


# MG sampler described by O'neill and Roberts in '98 paper
source("/home/macdona8/Documents/PhD/Simulating Epidemics/Metropolis-Gibbs (Oneill, Roberts).R")

# MG sampler described by Neal and Roberts in '03 paper
source("/home/macdona8/Documents/PhD/Simulating Epidemics/Metropolis-Gibbs (Neal, Roberts).R")


# ==== Test Runs ====

N = 20
beta = 0.25
gamma = 3

GSE.sim = GSE.simulation(N, a = 1, gamma, beta)
R = GSE.sim$removal.times
theta.gamma = 0.001
theta.beta = 0.001
theta = 0.001
no.its = 1000


MG.Neal = SIR.MetropolisGibbs(N, R, theta.gamma, theta.beta, no.its)

MG.Oneill = SIR.MetropolisGibbs2(N, R, theta.gamma, theta.beta, theta, no.its)

# ==== Experiments ====

# Vary Population Size

# Performance with increasing population 

# Varying infectious process and infectious period parameters (beta and gamma)
Sim = GSE.simulation(N = 20, a = 1, gamma = 3, beta = 0.25)

R = Sim$removal.times
I = Sim$infection.times
mean.Q1 = mean(R-I)
sum(R < 10000)
N = 20
theta.gamma = 0.001
theta.beta = 0.001
no.its = 3000