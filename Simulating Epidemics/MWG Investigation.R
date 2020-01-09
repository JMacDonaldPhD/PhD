#'
#' Metropolis-within-Gibbs Investigation (Neal-Roberts)
#' 
#'
#'
#'

# ==== Preamble ====

# Set Directory 

setwd("/home/macdona8/Documents/PhD/")
# MWG Algorithm
source("Simulating Epidemics/Metropolis-Gibbs (Neal, Roberts).R")


# Load in Simulated Datasets
source("Simulating Epidemics/Simulation for MWG investigation.R")

N = 200

low.finalsize = Low.Epidemic$final.size
low.remtimes = Low.Epidemic$removal.times

lowmed.finalsize = LowMed.Epidemic$final.size
lowmed.remtimes = LowMed.Epidemic$removal.times

med.finalsize = Med.Epidemic$final.size
med.remtimes =  Med.Epidemic$removal.times

highmed.finalsize = HighMed.Epidemic$final.size
highmed.remtimes = HighMed.Epidemic$removal.times
#simulation200 = readRDS("/home/macdona8/Documents/PhD/Simulating Epidemics/simulation.rds")
#simulation20 = readRDS("/home/macdona8/Documents/PhD/Simulating Epidemics/simulation20.rds")
#set.seed(2905)
#Bens.sim = GSE.simulation(N = 200, a = 1, gamma = 0.15, beta = 0.002, rinf.dist = "exp")

#N = Bens.sim$N
#R = Bens.sim$removal.times

#Low.Epi.MWG = SIR.MetropolisGibbs(N, R, gamma.prior = "exp", theta.gamma = 0.001,
#                                  beta.prior = "exp", theta.beta = 0.001, rinf.dist = "gamma", alpha = 1,
#                                  no.its = 10000, burn.in = 0, no.proposals = 1)

#plot(Low.Epi.MWG$draws[-(1:100),1], type = 'l')
#acf(Low.Epi.MWG$draws[-(1:100),1])


#plot(Low.Epi.MWG$draws[-(1:1000),1], type = 'l')
#acf(Low.Epi.MWG$draws[-(1:1000),1])
# ==== MWG Runs ====

# Low Final Size Epidemic (Final Size = 19)
no.proposals = round(i*low.finalsize)
Low.Epi.MWG = SIR.MetropolisGibbs(N = 200, low.remtimes, gamma.prior = "exp", theta.gamma = 0.001,
                               beta.prior = "exp", theta.beta = 0.001, rinf.dist = "exp",
                               no.its = 10000, burn.in = 0, no.proposals = 1)
par(mfrow = c(2,2), oma = c(0,0,2,0))
#title(main = "Epidemic with Low Final Size (19), N = 200, beta = 0.001, gamma = 0.15")
plot(Low.Epi.MWG$draws[-(1:1000),1], type = 'l', ylab = expression(beta))
plot(Low.Epi.MWG$draws[-(1:1000),2], type = 'l', ylab = expression(gamma))
acf(Low.Epi.MWG$draws[-(1:1000),1], main = expression(beta))
acf(Low.Epi.MWG$draws[-(1:1000),2], main = expression(gamma))
mtext('Epidemic with Low Final Size (19), N = 200, beta = 0.001, gamma = 0.15', outer = TRUE, cex = 1.5)


saveRDS(Low.Epi.MWG, paste(c("Centered MWG Investigation/Low.Epi.MWG_", no.proposals, ".rds"), sep = "", collapse = ""))

# Low-Medium Final Size (n_I = 73)

LowMed.Epi.MWG = SIR.MetropolisGibbs(N = 200, R = lowmed.remtimes, gamma.prior = "exp", theta.gamma = 0.001,
                                     beta.prior = "exp", theta.beta = 0.001, rinf.dist = "exp",
                                     no.its = 10000, burn.in = 0, no.proposals =1)
par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(LowMed.Epi.MWG$draws[-(1:1000),1], type = 'l', ylab = expression(beta))
plot(LowMed.Epi.MWG$draws[-(1:1000),2], type = 'l', ylab = expression(gamma))
acf(LowMed.Epi.MWG$draws[-(1:1000),1], main = expression(beta))
acf(LowMed.Epi.MWG$draws[-(1:1000),2], main = expression(gamma))
mtext('Epidemic with Low-Med Final Size (73), N = 200, beta = 0.001, gamma = 0.15', outer = TRUE, cex = 1.5)

saveRDS(LowMed.Epi.MWG, paste(c("Centered MWG Investigation/LowMed.Epi.MWG_", no.proposals, ".rds"), sep = "", collapse = ""))



# Medium Final Size (n_I = 101)
N = Med.Epidemic$N
R = Med.Epidemic$removal.times
Med.Epi.MWG = SIR.MetropolisGibbs(N, R, gamma.prior = "exp", theta.gamma = 0.001,
                                  beta.prior = "exp", theta.beta = 0.001, rinf.dist = "exp",
                                  no.its = 20000, burn.in = 0, no.proposals = 1)

par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(Med.Epi.MWG$draws[-(1:1000),1], type = 'l', ylab = expression(beta))
plot(Med.Epi.MWG$draws[-(1:1000),2], type = 'l', ylab = expression(gamma))
acf(Med.Epi.MWG$draws[-(1:1000),1], main = expression(beta))
acf(Med.Epi.MWG$draws[-(1:1000),2], main = expression(gamma))
mtext('Epidemic with Medium Final Size (101), N = 200, beta = 0.001, gamma = 0.15', outer = TRUE, cex = 1.5)

saveRDS(Med.Epi.MWG, paste(c("Centered MWG Investigation/Med.Epi.MWG_", no.proposals, ".rds"), sep = "", collapse = ""))


# High-Medium Final Size (n_I = 125)
N = HighMed.Epidemic$N
R = HighMed.Epidemic$removal.times
HighMed.Epi.MWG = SIR.MetropolisGibbs(N, highmed.remtimes, gamma.prior = "exp", theta.gamma = 0.001,
                                      beta.prior = "exp", theta.beta = 0.001, rinf.dist = "exp",
                                      no.its = 20000, burn.in = 0, no.proposals = 10)


par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(HighMed.Epi.MWG$draws[-(1:1000),1], type = 'l', ylab = expression(beta))
plot(HighMed.Epi.MWG$draws[-(1:1000),2], type = 'l', ylab = expression(gamma))
acf(HighMed.Epi.MWG$draws[-(1:1000),1], main = expression(beta))
acf(HighMed.Epi.MWG$draws[-(1:1000),2], main = expression(gamma))
mtext('Epidemic with High-Med Final Size (125), N = 200, beta = 0.001, gamma = 0.15', outer = TRUE, cex = 1.5)

saveRDS(HighMed.Epi.MWG, paste(c("Centered MWG Investigation/HighMed.Epi.MWG_", no.proposals, ".rds"), sep = "", collapse = ""))

saveRDS(test.run, "test.propose20.1.rds")









