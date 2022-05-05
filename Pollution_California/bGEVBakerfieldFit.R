##############################################################################
## bGEV model fit for the Bakerfield NO2 data, selected using WAIC/DIC/PITs ##
##############################################################################
## Data and libraries
library(INLA)
source('Pollution_California/utilsBakersfield.R')
inla.data = read.table('Pollution_California/Bakersfield.txt', header = T)
dim(inla.data)
sum(!is.na(inla.data$NO2_max))
names(inla.data)
head(inla.data)
length(unique(inla.data$year))
inla.data$year = inla.data$year - 1999 
inla.data$month = inla.group(inla.data$month, method="cut",n = 12)


## Setup
a           = 0.01
u           = sd(inla.data$NO2_max, na.rm = T)
cyclic      = T
hyperRW     = list(theta = list(prior = "pc.prec", param = c(u, a)))
hyper.month = list(prec = list(initial = log(1/0.024), fixed = TRUE))

## No model for tail and spread
null.matrix = matrix(nrow = nrow(inla.data), ncol = 0)
spread.x    = null.matrix
tail.x      = null.matrix

## INLA setup
control.family    = list(hyper = hyper.bgev, control.bgev = control.bgev)
control.predictor = list(compute = TRUE)
control.compute   = list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)
control.inla      = list(strategy = "simplified.laplace", int.strategy = "eb")

## Model fit
formula = inla.mdata(NO2_max, spread.x, tail.x) ~ -1 + intercept + year + f(month, model = 'rw2', hyper = hyper.month, cyclic = cyclic, scale.model = TRUE) + WindS_max +
  f(Temp_meanL, model = 'rw2', hyper = hyperRW, cyclic = cyclic, scale.model = TRUE)
fit10    = inla(formula,
                family = "bgev",
                data   = inla.data,
                control.family    = control.family,
                control.predictor = control.predictor,
                control.compute   = control.compute,
                control.inla      = control.inla,
                verbose           = F)
