#'
#' MWG Debug (*NEW*)
#'
#'
#'

# ==== Preamble ====

source("/home/macdona8/Documents/PhD/Simulating Epidemics/Metropolis-Gibbs (Neal, Roberts).R")

# ==== Simulating Test Data ====
# Set Seed to ensure we get an epidemic which takes off
set.seed(1)


Sim = GSE.simulation(N = 200, a = 1, gamma = 0.15, beta = 0.001, rinf.dist = "exp")

# Extract output from simulation for use in MWG sampler
Sim$final.size
N = Sim$N
R = Sim$removal.times
I = Sim$infection.times

# ==== Running Algortihm ====

Test = SIR.MetropolisGibbs(N, R, theta.gamma = 0.001, theta.beta = 0.001, beta.prior = "exp",
                    gamma.prior = "exp", rinf.dist = "exp", alpha = 1, no.its = 10000, 
                    burn.in = 0, no.proposals = 1)

