#############################################################################
## Negative log-likelihood functions associates to the GEV and bGEV models ##
#############################################################################
# By Daniela C.C.

source('Code/bGEVcode.R')
source('Code/GEVcode.R')


# Neg log-lik of GEV using classical parametrisation
nllik.gev = function(par, x, log = TRUE){
  mu    = par[1]
  sigma = par[2]
  xi    = par[3]
  ll    = rep(NA, length(x))
  for(i in 1:length(x))
    ll[i] = dgev(x[i], mu, sigma, xi, log = log)
  -sum(ll)
}

# Neg log-lik of bGEV using classical parametrisation
nllik.bgev = function(par, x, p_a, p_b, s, log = TRUE){
  mu    = par[1]
  sigma = par[2]
  xi    = par[3]
  ll    = rep(NA, length(x))
  if(xi < 0){
    return(1e10)
  }else{
    for(i in 1:length(x))
      ll[i] = dbgev(x[i], mu, sigma, xi, p_a = p_a, p_b = p_b, s = s, log = log)
    return(-sum(ll))
  }
}

# Neg log-lik of GEV using classical parametrisation with covariate-dependent location
nllik.gevx = function(par, x, w, log = TRUE){
  mu0   = par[1]
  mu1   = par[2]
  sigma = par[3]
  xi    = par[4]
  mu    = mu0 + mu1*w
  ll    = rep(NA, length(x))
  for(i in 1:length(x))
    ll[i] = dgev(x[i], mu[i], sigma, xi, log = log)
    
  -sum(ll)
}

# Neg log-lik of bGEV using classical parametrisation with covariate-dependent location
nllik.bgevx = function(par, x, w, p_a, p_b, s, log = TRUE){
  mu0   = par[1]
  mu1   = par[2]
  sigma = par[3]
  xi    = par[4]
  mu    = mu0 + mu1*w
  ll    = rep(NA, length(x))
  if(xi < 0){
    return(1e10)
  }else{
    for(i in 1:length(x))
      ll[i] = dbgev(x[i], mu[i], sigma, xi, p_a = p_a, p_b = p_b, s = s, log = log)
    return(-sum(ll))
  }
}

# Neg log-lik of GEV using new parametrisation 
nllik.gev2 = function(par, x, alpha = 0.5, beta = 0.5, log = TRUE){
  q  = par[1]
  s  = par[2]
  xi = par[3]
  ll = rep(NA, length(x))
  if(xi < 0){
    return(1e10)
  }else{
    for(i in 1:length(x))
      ll[i] = dgev2(x[i], q = q, s = s, xi = xi, alpha = alpha, beta = beta, log = log)
    return(-sum(ll))
  }
}

# Neg   log-lik of bGEV using new parametrisation
nllik.bgev2 = function(par, x, alpha = 0.5, beta = 0.5, p_a, p_b, s, log = TRUE){
  tmp   = new.to.old(par, alpha = alpha, beta = beta)
  mu    = tmp$mu
  sigma = tmp$sigma 
  ll    = rep(NA, length(x))
  if(xi < 0){
    return(1e10)
  }else{
    for(i in 1:length(x))
      ll[i] = dbgev(x[i], mu, sigma, xi, p_a = p_a, p_b = p_b, s = s, log = log)
    return(-sum(ll))
  }
}
