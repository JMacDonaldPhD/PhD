#'
#' Simulating a Poisson Process
#' 

#'
#' Simulate the number of events observed in
#' a homogenous poisson process in 
#' period (0,t) by drawing a poisson random 
#' variable with rate lambda*t
#'
lambda = 1
t = 10
n = rpois(1, lambda = lambda*t)

#'
#' As shown in the 'Epidemic Notes'
#' the arrival times in a poisson process are 
#' distributed as N(t) = n iid uniform random variables
#' on the interval (0,t)

S = sort(runif(n, min = 0, max = t))

#'
#' Bringing this all into a function,
#'

pois.process = function(rate, t){
  n = rpois(1, lambda = rate*t)
  S = sort(runif(n, min = 0, max = t))
  return(list(events = n, times = S))
}

pois.process(1, 10)
