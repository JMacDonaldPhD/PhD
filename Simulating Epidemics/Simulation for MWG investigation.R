#'
#' Simulating An Epidemic for use with MWG Algorithm
#'
#'

# Simulate Epidemics with a closed population of N = 200 individuals

# Want Epidemics of various final size, i.e how many individuals are infected
# as a result of one individual being infected at the start.

# Low, low-med, high-med, high final size


# ==== Small Scale Epidemic vs. Large Scale Epidemic ====
# Size 20

final.size = 0
while(final.size == 0){
  LowPopn.Epidemic = GSE.simulation(N = 20, a = 1, gamma = 0.15, beta = 0.01, rinf.dist = "exp")
  final.size = LowPopn.Epidemic$final.size
}

# Size 200

final.size = 0
while(final.size == 0){
  MedPopn.Epidemic = GSE.simulation(N = 200, a  = 1, gamma = 0.15, beta = 0.001, rinf.dist = "exp")
  final.size = MedPopn.Epidemic$final.size
}

# Size 2000

final.size = 0
while(final.size == 0){
  HighPopn.Epidemic = GSE.simulation(N = 2000, a  = 1, gamma = 0.15, beta = 0.1, rinf.dist = "exp")
  final.size = HighPopn.Epidemic$final.size
}

# Size 200
# set the seed 100 == Low
#              3   == Medium 
#              11  == Low-Medium 
#              93  == High-Medium

#set.seed(100)
#Low.Epidemic = GSE.simulation(N = 200, a  = 1, gamma = 0.15, beta = 0.001, rinf.dist = "exp")

#set.seed(3)
#Med.Epidemic = GSE.simulation(N = 200, a  = 1, gamma = 0.15, beta = 0.001, rinf.dist = "exp")

#set.seed(11)
#LowMed.Epidemic = GSE.simulation(N = 200, a  = 1, gamma = 0.15, beta = 0.001, rinf.dist = "exp")

#set.seed(93)
#HighMed.Epidemic = GSE.simulation(N = 200, a  = 1, gamma = 0.15, beta = 0.001, rinf.dist = "exp")


#low = c(0,0)
#high =  c(0,0)
#highmed =  c(0,0)
#for(i in 1:500){
#  set.seed(i)
#  simulation200 = GSE.simulation(N = 200, a  = 1, gamma = 0.15, beta = 0.001, rinf.dist = "exp")
#  FS = simulation200$final.size 
#  if(175 < FS){
#    high = cbind(high, c(i, FS))
#  } else if(FS > 0 & FS < 25){
#    low = cbind(low, c(i, FS))
#  } else if(FS < 150 & FS > 110){
#    highmed = cbind(highmed, c(i, FS))
#  }
#}

# ==== Exporting Simulations ====

#saveRDS(simulation200, "/home/macdona8/Documents/PhD/Simulating Epidemics/simulation200.rds")
#saveRDS(simulation20, "/home/macdona8/Documents/PhD/Simulating Epidemics/simulation20.rds")
