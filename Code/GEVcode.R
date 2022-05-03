#################################################
## Functions to work with the GEV distribution ##
#################################################
# By Silius M.V. and Daniela C.C.
source('Code/utils.R')

# PDF
pgev = function(x, mu, sigma, xi) {
  fix_lengths(x, mu, sigma, xi)
  ifelse(xi == 0,
         exp(-exp(- (x - mu) / sigma)),
         exp(-pmax(0, 1 + xi * (x - mu) / sigma) ^ (-1 / xi)))
}

# Quantile function
qgev = function(p, mu, sigma, xi) {
  fix_lengths(p, mu, sigma, xi)
  ifelse(xi == 0,
         mu - sigma * log(-log(p)),
         mu - sigma * (1 / xi) * (1 - (- log(p)) ^ (-xi)))
}

# Random generation
rgev = function(n, mu, sigma, xi) {
  lengths = sapply(list(mu, sigma, xi), length)
  if (any(lengths > n)) stop("Bad input lengths")
  qgev(runif(n), mu, sigma, xi)
}

# Density using the usual parametrisation
dgev = function(x, mu, sigma, xi, log = FALSE) {
  fix_lengths(x, mu, sigma, xi)
  res = ifelse(xi == 0,
               -exp(- (x - mu) / sigma),
               -pmax(0, 1 + xi * (x - mu) / sigma) ^ (-1 / xi))
  res = res - log(sigma) +
    ifelse(xi == 0,
           - (x - mu) / sigma,
           ifelse(1 + xi * (x - mu) / sigma > 0,
                  - (1 / xi + 1) * log(1 + xi * (x - mu) / sigma),
                  -Inf))
  if (!log) res = exp(res)
  res
}

# Density using the new parametrisation
dgev2 = function(x, q, sb, xi, alpha = 0.5, beta = 0.5, log = FALSE) {
  tmp   = new.to.old(c(q,sb,xi), alpha = alpha, beta = beta)
  mu    = tmp$mu
  sigma = tmp$sigma 
  
  res = ifelse(xi == 0,
               -exp(- (x - mu) / sigma),
               -pmax(0, 1 + xi * (x - mu) / sigma) ^ (-1 / xi))
  res = res - log(sigma) +
    ifelse(xi == 0,
           - (x - mu) / sigma,
           ifelse(1 + xi * (x - mu) / sigma > 0,
                  - (1 / xi + 1) * log(1 + xi * (x - mu) / sigma),
                  -Inf))
  if (!log) res = exp(res)
  res
}

# Return levels
return_level_gev = function(period, mu, sigma, xi) {
  if (any(period <= 1)) warning("invalid period")
  p = ifelse(period > 1, 1 - 1 / period, NA)
  qgev(p, mu, sigma, xi)
}
