############################################################
##                                                        ##
## Sim1: compare bGEV and GEV at estimating return levels ##
##                                                        ##
############################################################
source('Code/GEVbGEVlikelihoods.R')
library(rje)
library(xtable)
## Setup
M     = 500
mu    = 0
sigma = 1
xi    = 0.1
alpha = 1/xi
p_a   = 0.05
p_b   = 0.2
s     = 5
rp    = c(30,50,100)
rl    = evd::qgev(1/rp, mu, sigma, xi, lower.tail = F) # true rls

## Simulation scenarios
n = c(30,50,100,500)
N = c(30,50,100,500,1000)
scenarios = expand.grid(n, N); colnames(scenarios) = c('n', 'N'); scenarios

rl.bgev30 = rl.bgev50 = rl.bgev100 = NULL
rl.gev30  = rl.gev50 = rl.gev100 = NULL
t0 = Sys.time()
for(mc in 5:M){
  printPercentage(mc, M)
  tmp.bgev30 = tmp.bgev50 = tmp.bgev100 = tmp.gev30 = tmp.gev50 = tmp.gev100 = NULL
  
  for(i in 1:nrow(scenarios)){
    print(paste('Scenario', i))
    N = scenarios$N[i]
    n = scenarios$n[i]
    y = rep(NA, N)
    for(j in 1:N){
      set.seed(n*j+mc)
      y[j] = max(evd::rfrechet(n, shape = alpha))
    }
    in2  = sqrt(6 * var(y))/pi
    in1  = mean(y) - 0.57722 * in2
    init = c(in1, in2, 0.1)
    
    fit.bgev = tryCatch(optim(init, nllik.bgev, x = y, p_a = p_a, p_b = p_b, s = s), error = function(e) e)
    fit.gev  = tryCatch(optim(init, nllik.gev, x = y), error = function(e) e)
    
    if(!inherits(fit.bgev, 'error') & !inherits(fit.gev, 'error')){
      tmp1 = return_level_bgev(period = rp, 
                               mu = fit.bgev$par[1], sigma = fit.bgev$par[2], xi = fit.bgev$par[3], 
                               p_a = p_a, p_b = p_b, s = s)
      tmp2 = return_level_bgev(period = rp, 
                               mu = fit.gev$par[1], sigma = fit.gev$par[2], xi = fit.gev$par[3])
      
      tmp.bgev30  = c(tmp.bgev30, tmp1[1])
      tmp.bgev50  = c(tmp.bgev50, tmp1[2])
      tmp.bgev100 = c(tmp.bgev100, tmp1[3])
      
      tmp.gev30   = c(tmp.gev30, tmp2[1])
      tmp.gev50   = c(tmp.gev50, tmp2[2])
      tmp.gev100  = c(tmp.gev100, tmp2[3])
      
      
    }
  }
  rl.bgev30  = rbind(rl.bgev30, tmp.bgev30)
  rl.bgev50  = rbind(rl.bgev50, tmp.bgev50)
  rl.bgev100 = rbind(rl.bgev100, tmp.bgev100)
  
  rl.gev30   = rbind(rl.gev30, tmp.gev30)
  rl.gev50   = rbind(rl.gev50, tmp.gev50)
  rl.gev100  = rbind(rl.gev100, tmp.gev100)
  
  write.table(rl.bgev30, file = 'Simulation/Sim1/sim1out_bgev_rl30.txt', sep = '\t', row.names = F, col.names = F)
  write.table(rl.bgev50, file = 'Simulation/Sim1/sim1out_bgev_rl50.txt', sep = '\t', row.names = F, col.names = F)
  write.table(rl.bgev100, file = 'Simulation/Sim1/sim1out_bgev_rl100.txt', sep = '\t', row.names = F, col.names = F)
  
  write.table(rl.gev30, file = 'Simulation/Sim1/sim1out_gev_rl30.txt', sep = '\t', row.names = F, col.names = F)
  write.table(rl.gev50, file = 'Simulation/Sim1/sim1out_gev_rl50.txt', sep = '\t', row.names = F, col.names = F)
  write.table(rl.gev100, file = 'Simulation/Sim1/sim1out_gev_rl100.txt', sep = '\t', row.names = F, col.names = F)
}
Sys.time() - t0
# hist(y, freq = F); curve(evd::dfrechet(x, mu, n^(1/alpha)*sigma, alpha), from = 0, to = 10, add = T)

## RMSE bgev
rl.bgev30  = read.table('Simulation/Sim1/sim1out_bgev_rl30.txt', header = F); dim(rl.bgev30)
rl.bgev50  = read.table('Simulation/Sim1/sim1out_bgev_rl50.txt', header = F); dim(rl.bgev50)
rl.bgev100 = read.table('Simulation/Sim1/sim1out_bgev_rl100.txt', header = F); dim(rl.bgev100)

rmse.bgev30  = round(sqrt(colMeans((rl.bgev30 - rl[1])^2)),2)
rmse.bgev50  = round(sqrt(colMeans((rl.bgev50 - rl[2])^2)),2)
rmse.bgev100 = round(sqrt(colMeans((rl.bgev100 - rl[3])^2)),2)

rmse.bgev30           = data.frame(matrix(rmse.bgev30, ncol = length(N)))
rmse.bgev50           = data.frame(matrix(rmse.bgev50, ncol = length(N)))
rmse.bgev100           = data.frame(matrix(rmse.bgev100, ncol = length(N)))
colnames(rmse.bgev30) = colnames(rmse.bgev50) = colnames(rmse.bgev100) = N
rownames(rmse.bgev30) = rownames(rmse.bgev50) = rownames(rmse.bgev100) = n
rmse.bgev30
rmse.bgev50
rmse.bgev100

## RMSE gev
rl.gev30  = read.table('Simulation/Sim1/sim1out_gev_rl30.txt', header = F); dim(rl.gev30)
rl.gev50  = read.table('Simulation/Sim1/sim1out_gev_rl50.txt', header = F); dim(rl.gev50)
rl.gev100 = read.table('Simulation/Sim1/sim1out_gev_rl100.txt', header = F); dim(rl.gev100)

rmse.gev30  = round(sqrt(colMeans((rl.gev30 - rl[1])^2)),2)
rmse.gev50  = round(sqrt(colMeans((rl.gev50 - rl[2])^2)),2)
rmse.gev100 = round(sqrt(colMeans((rl.gev100 - rl[3])^2)),2)

rmse.gev30           = data.frame(matrix(rmse.gev30, ncol = length(N)))
rmse.gev50           = data.frame(matrix(rmse.gev50, ncol = length(N)))
rmse.gev100          = data.frame(matrix(rmse.gev100, ncol = length(N)))
colnames(rmse.gev30) = colnames(rmse.gev50) = colnames(rmse.gev100) = N
rownames(rmse.gev30) = rownames(rmse.gev50) = rownames(rmse.gev100) = n
rmse.gev30
rmse.gev50
rmse.gev100

## Difference
rmse30  = rmse.gev30 -rmse.bgev30; rmse30
rmse50  = rmse.gev50 -rmse.bgev50; rmse50
rmse100 = rmse.gev100 -rmse.bgev100; rmse100

rmse        = data.frame(row.names = rownames(rmse.gev30))
rmse$`30`   = paste0(round(rmse30[,1],2), '/', round(rmse50[,1],2), '/', round(rmse100[,1],2))
rmse$`50`   = paste0(round(rmse30[,2],2), '/', round(rmse50[,2],2), '/', round(rmse100[,2],2))
rmse$`100`  = paste0(round(rmse30[,3],2), '/', round(rmse50[,3],2), '/', round(rmse100[,3],2))
rmse$`500`  = paste0(round(rmse30[,4],2), '/', round(rmse50[,4],2), '/', round(rmse100[,4],2))
rmse$`1000` = paste0(round(rmse30[,5],2), '/', round(rmse50[,5],2), '/', round(rmse100[,5],2))
xtable(rmse)










