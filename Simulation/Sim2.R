#################################################################################################
##                                                                                             ##
## Sim2: assess the effect of the hyperparameters $p_a, p_b$ and $c_1=c_2$ over the bGEV model ##
##                                                                                             ##
#################################################################################################
source('Code/GEVbGEVlikelihoods.R')
library(rje)
## Setup
N     = 100
M     = 500
mu    = 0
sigma = 1
xi    = 0.2
init  = c(mu, sigma, xi)
rp    = 50
rl    = evd::qgev(1/rp, mu, sigma, xi, lower.tail = F) # true 50-year rl

## Simulation scenarios
p_a       = c(0.05, 0.1, 0.15)
p_b       = c(0.2, 0.25, 0.3)
s         = c(3, 5)
scenarios = expand.grid(p_a, p_b, s); colnames(scenarios) = c('pa', 'pb', 's'); scenarios

## MC replicates
rl.bgev = rl.gev = NULL
t0 = Sys.time()
for(mc in 1:M){
  printPercentage(mc, M)
  
  tmp.bgev = tmp.gev = NULL
  set.seed(N*mc)
  y = rgev(N, mu, sigma, xi)
  
  for(i in 1:nrow(scenarios)){
    print(paste('Scenario', i))
    
    hyperpar = scenarios[i, ]
    fit.bgev = tryCatch(optim(init, nllik.bgev, x = y, p_a = hyperpar$pa, p_b = hyperpar$pb, s = hyperpar$s), error = function(e) e)
    if(!inherits(fit.bgev, 'error')){
      tmp.bgev = c(tmp.bgev, return_level_bgev(period = rp, 
                                               mu = fit.bgev$par[1], sigma = fit.bgev$par[2], xi = fit.bgev$par[3], 
                                               p_a = hyperpar$pa, p_b = hyperpar$pb, s = hyperpar$s))
      
    }
  }
  fit.gev = tryCatch(optim(init, nllik.gev, x = y), error = function(e) e)
  if(!inherits(fit.bgev, 'error')){
    tmp.gev = return_level_gev(period = rp, 
                               mu = fit.gev$par[1], sigma = fit.gev$par[2], xi = fit.gev$par[3])
  }
  rl.bgev = rbind(rl.bgev, tmp.bgev)
  rl.gev = rbind(rl.gev, tmp.gev)
  write.table(rl.bgev, file = 'Simulation/Sim2/sim2out_bgev.txt', sep = '\t', row.names = F, col.names = F)
  write.table(rl.gev, file = 'Simulation/Sim2/sim2out_gev.txt', sep = '\t', row.names = F, col.names = F)
}
Sys.time() - t0 #  5.4 hours

## RMSE bgev
rl.bgev   = read.table('Simulation/Sim2/sim2out_bgev.txt', header = F)
dim(rl.bgev)
rmse.bgev = round(sqrt(colMeans((rl.bgev - rl)^2)),2)
tmp = cbind(scenarios, c(rmse.bgev))

rmse.bgev = matrix(c(tmp[1,4],  tmp[4,4],  tmp[7,4],
                     tmp[10,4], tmp[13,4], tmp[16,4],
                     tmp[2,4],  tmp[5,4],  tmp[8,4],
                     tmp[11,4], tmp[14,4], tmp[17,4],
                     tmp[3,4],  tmp[6,4], tmp[9,4],
                     tmp[12,4], tmp[15,4], tmp[18,4]),
                   ncol = length(p_a)*length(s))
rownames(rmse.bgev) = p_b
xtable(rmse.bgev)

## RMSE gev
rl.gev   = read.table('Simulation/Sim2/sim2out_gev.txt', header = F)$V1
rmse.gev = round(sqrt(mean((rl.gev - rl)^2)),2); rmse.gev
