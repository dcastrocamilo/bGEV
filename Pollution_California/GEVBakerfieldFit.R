#############################################################
## GEV fit for the bGEV model selected using WAIC/DIC/PITs ##
#############################################################

library(INLA)
source('Pollution_California/utilsBakersfield.R')
source('Code/GEVcode.R')
source('Code/bGEVcode.R')

## Data
inla.data = read.table('Pollution_California/Bakersfield.txt', header = T)
inla.data$year = inla.data$year - 1999 
inla.data$month = inla.group(inla.data$month, method="cut",n = 12)


## Setup
a           = 0.01
u           = sd(inla.data$NO2_max, na.rm = T)
s           = 0.1 # fixed scaling for GEV precision
sxi         = 0.01 # fixed scaling for GEV shape (gev.scale.xi)
rp          = c(2,30, 50, 100)*12
cyclic      = T
hyperRW     = list(theta = list(prior = "pc.prec", param = c(u, a)))
hyper.month = list(prec = list(initial = log(1/0.024), fixed = TRUE))

## INLA setup
hyper.prec        = list(initial = 1, fixed = F, prior = "loggamma", param = c(3, 3))
hyper.tail        = list(initial = 0.5, fixed = F, prior = "loggamma", param = c(1, 9)) # chosen to be close to pc.gevtail(lambda=7)
hyper.gev         = list(theta1 = hyper.prec, theta2 = hyper.tail)
control.family    = list(hyper = hyper.gev, gev.scale.xi = sxi)
control.predictor = list(compute = TRUE)
control.compute   = list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)
control.inla      = list(strategy = "simplified.laplace", int.strategy = "eb")

## Model formula
formula = NO2_max ~ -1 + intercept + year + f(month, model = 'rw2', hyper = hyper.month, cyclic = cyclic, scale.model = TRUE) + WindS_max +
  f(Temp_meanL, model = 'rw2', hyper = hyperRW, cyclic = cyclic, scale.model = TRUE)


## GEV fit
# bgev.start = fit.bgev$mode$theta
# bgev.start[2] = map.tail(bgev.start[2], interval = tail.interval, inverse = F)/sxi
fit.gev = inla(formula,
               family = "gev",
               data   = inla.data,
               # control.mode = list(theta = bgev.start, restart = TRUE), 
               control.family    = control.family,
               control.predictor = control.predictor,
               control.compute   = control.compute,
               control.inla      = control.inla,
               verbose           = T) 

# save(fit.gev, file = 'Pollution_California/gevfit.Rdata')
