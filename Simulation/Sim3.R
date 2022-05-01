#############################################################################################
##                                                                                         ##
## Sim3: assess the effect of the hyperparameters $\alpha$ and $\beta$ over the bGEV model ##
##                                                                                         ##
#############################################################################################
source('Code/GEVbGEVlikelihoods.R')
library(rje)
library(xtable)
## Setup
N     = 100
M     = 500
mu    = 0
sigma = 1
xi    = 0.2
init  = c(mu, sigma, xi)
p_a   = 0.05
p_b   = 0.2
s     = 5
rp    = 50
rl    = evd::qgev(1/rp, mu, sigma, xi, lower.tail = F) # true 50-year rl

## Simulation scenarios
alpha     = c(0.1, 0.3, 0.5, 0.7, 0.9) 
beta      = c(0.05, 0.3, 0.5, 0.7, 0.9)
scenarios = expand.grid(alpha, beta); colnames(scenarios) = c('alpha', 'beta'); scenarios

## MC replicates
rl.bgev = rl.gev = NULL
t0 = Sys.time()
for(mc in 1:M){
  printPercentage(mc, M)
  tmp.bgev = NULL
  set.seed(N*mc)
  y = rgev(N, mu, sigma, xi)
  
  for(i in 1:nrow(scenarios)){
    print(paste('Scenario', i))
    
    hyperpar = scenarios[i, ]
    init2    = as.numeric(old.to.new(init, alpha = hyperpar$alpha, beta = hyperpar$beta))
    fit.bgev = tryCatch(optim(init2, nllik.bgev2, x = y, alpha = hyperpar$alpha, beta = hyperpar$beta, p_a = p_a, p_b = p_b, s = s), error = function(e) e)
    
    if(!inherits(fit.bgev, 'error')){
      par      = as.numeric(new.to.old(fit.bgev$par, alpha = hyperpar$alpha, beta = hyperpar$beta))
      tmp.bgev = c(tmp.bgev, return_level_bgev(period = rp, 
                                               mu = par[1], sigma = par[2], xi = par[3], 
                                               p_a = p_a, p_b = p_a, s = s))
    }
  }
  fit.gev  = tryCatch(optim(init, nllik.gev, x = y), error = function(e) e)
  if(!inherits(fit.bgev, 'error')){
    tmp.gev = return_level_gev(period = rp, 
                               mu = fit.gev$par[1], sigma = fit.gev$par[2], xi = fit.gev$par[3])
    rl.gev  = rbind(rl.gev, tmp.gev)
  }
  rl.bgev = rbind(rl.bgev, tmp.bgev)
  write.table(rl.bgev, file = 'Simulation/Sim3/sim3out_bgev.txt', sep = '\t', row.names = F, col.names = F)
  write.table(rl.gev, file = 'Simulation/Sim3/sim3out_gev.txt', sep = '\t', row.names = F, col.names = F)
}
Sys.time() - t0

## RMSE bgev
rl.bgev             = read.table('Simulation/Sim3/sim3out_bgev.txt', header = F); dim(rl.bgev)
rmse.bgev           = round(sqrt(colMeans((rl.bgev - rl)^2)),2); cbind(scenarios, c(rmse.bgev))
rmse.bgev           = data.frame(matrix(rmse.bgev, ncol = length(beta)))
colnames(rmse.bgev) = beta
rownames(rmse.bgev) = alpha
rmse.bgev = rmse.bgev[-1,-(1:2)]
xtable(rmse.bgev)

## RMSE gev
rl.gev   = read.table('Simulation/Sim3/sim3out_gev.txt', header = F)$V1
rmse.gev = round(sqrt(mean((rl.gev - rl)^2)),2); rmse.gev








