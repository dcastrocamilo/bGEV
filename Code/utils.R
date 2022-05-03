#################################################
## Utility functions to work with the bGEV-GEV ##
#################################################
# By Silius M.V. and Daniela C.C.

new.to.old = function(par, alpha = 0.5, beta = 0.5){
  q = par[1] 
  s = par[2]
  xi = par[3]
  if(xi == 0){
    ell1 = log(-log(alpha))
    ell2 = log(-log(beta/2))
    ell3 = log(-log(1-beta/2))
  }else{
    ell1 = (-log(alpha))^(-xi)
    ell2 = (-log(beta/2))^(-xi)
    ell3 = (-log(1-beta/2))^(-xi)
  }
  if(xi == 0){
    sigma = s/(ell2 - ell3)
    mu    = q + sigma*ell1
  }else{
    mu    = q-s*(ell1 - 1)/(ell3-ell2)
    sigma = xi*s/(ell3-ell2)
  }
  list(mu = mu, sigma = sigma, xi = xi)
}
old.to.new = function(par, alpha = 0.5, beta = 0.5){
  mu    = par[1]
  sigma = par[2]
  xi    = par[3]
  
  qalpha = qgev(alpha, mu, sigma, xi)
  qbeta1 = qgev(beta/2, mu, sigma, xi)
  qbeta2 = qgev((1-beta/2), mu, sigma, xi)
  list(q = qalpha, s = qbeta2 - qbeta1, xi = xi)
}

fix_lengths = function(...) {
  call = match.call()
  varnames = sapply(call[-1], as.character)
  e = parent.frame()
  vars = lapply(varnames, get, envir = e)
  lengths = sapply(vars, length)
  max_length = max(lengths)
  if (any(max_length %% lengths != 0)) stop("Bad input lengths")
  for (i in seq_along(vars)) {
    if (lengths[i] < max_length) {
      assign(varnames[i], rep(vars[[i]], max_length / lengths[i]), envir = e)
    }
  }
  0
}