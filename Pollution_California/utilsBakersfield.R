################################################################################
## Utility functions for the INLA bGEV model applied to NO2 at Bakerfield, CA ##
################################################################################

## Utility function to define hyper for tail
map.tail = function(x, interval, inverse = FALSE) {
  if (!inverse) {
    return (interval[1] + (interval[2] - interval[1]) * exp(x)/(1.0 + exp(x)))
  } else {
    return (log((x-interval[1])/(interval[2]-x)))
  }
}

## Default hyper for spread 
hyper.spread = list(initial = 1,
                    fixed   = FALSE,
                    prior   = "loggamma",
                    param   = c(3, 3))

## Default hyper for tail 
tail.interval = c(0, 0.5)
lambda        = 7
tail          = 0.08
tail.intern   = map.tail(tail, tail.interval, inverse=TRUE)
hyper.tail    = list(initial = if (tail == 0.0) -Inf else tail.intern, 
                     prior = "pc.gevtail",
                     param = c(lambda, tail.interval), 
                     fixed = if (tail == 0.0) TRUE else FALSE)

## All hypers
hyper.bgev = list(spread = hyper.spread,
                  tail   = hyper.tail)

## bGEV hyperparameters 
alpha = 0.5  # alpha for q_alpha
beta  = 0.25 # beta for s_beta
pa    = 0.05 # Probability pa for the mixing region
pb    = 0.2  # Probability pb for the mixing region
c1    = 5    # parameters of the Beta weight function. Hard-coded for the moment: c1=c2=5

## INLA control
control.bgev = list(q.location = alpha,
                    q.spread   = beta,
                    q.mix      = c(pa, pb), 
                    beta.ab    = c1)
