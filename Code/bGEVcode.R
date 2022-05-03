##################################################
## Functions to work with the bGEV distribution ##
##################################################
# By Silius M.V. and Daniela C.C.
source('Code/utils.R')

# PDF
pbgev = function(x, mu, sigma, xi, p_a = .1, p_b = .2, s = 5) {
  fix_lengths(x, mu, sigma, xi, p_a, p_b, s)
  g = get_gumbel_par(mu, sigma, xi, p_a, p_b)
  a = qgev(p_a, mu, sigma, xi)
  b = qgev(p_b, mu, sigma, xi)
  p = pbeta((x - a) / (b - a), s, s)
  pgev(x, mu, sigma, xi) ^ p * pgev(x, g$mu, g$sigma, 0) ^ (1 - p)
}

# Quantile function
qbgev = function(p, mu, sigma, xi, p_a = .1, p_b = .2, s = 5) {
  fix_lengths(p, mu, sigma, xi, p_a, p_b, s)
  res = rep(NA, length(p))
  gumbel = which(p <= p_a)
  frechet = which(p >= p_b)
  mixing = which(p_a < p & p < p_b)
  if (any(gumbel)) {
    g = get_gumbel_par(mu[gumbel], sigma[gumbel], xi[gumbel], p_a[gumbel], p_b[gumbel])
    res[gumbel] = qgev(p[gumbel], g$mu, g$sigma, 0)
  }
  if (any(frechet)) res[frechet] = qgev(p[frechet], mu[frechet], sigma[frechet], xi[frechet])
  if (any(mixing)) {
    res[mixing] = qbgev_mixing(p[mixing], mu[mixing], sigma[mixing], xi[mixing],
                               p_a[mixing], p_b[mixing], s[mixing])
  }
  res
}

# Random generation
rbgev = function(n, mu, sigma, xi, p_a = .1, p_b = .2, s = 5) {
  lengths = sapply(list(mu, sigma, xi), length)
  if (any(lengths > n)) stop("Bad input lengths")
  qbgev(runif(n), mu, sigma, xi, p_a, p_b, s)
}

# Density using the usual parametrisation
dbgev = function(x, mu, sigma, xi, p_a = .1, p_b = .2, s = 5, log = FALSE) {
  fix_lengths(x, mu, sigma, xi, p_a, p_b, s)
  a = qgev(p_a, mu, sigma, xi)
  b = qgev(p_b, mu, sigma, xi)
  res = rep(NA, length(x))
  gumbel = which(x <= a)
  frechet = which(x >= b)
  mixing = which(a < x & x < b)
  if (any(gumbel)) {
    g = get_gumbel_par(mu[gumbel], sigma[gumbel], xi[gumbel], p_a[gumbel], p_b[gumbel])
    res[gumbel] = dgev(x[gumbel], g$mu, g$sigma, 0, log = log)
  }
  if (any(frechet)) res[frechet] = dgev(x[frechet], mu[frechet], sigma[frechet], xi[frechet], log = log)
  if (any(mixing)) {
    res[mixing] = dbgev_mixing(x[mixing], mu[mixing], sigma[mixing], xi[mixing],
                               p_a[mixing], p_b[mixing], s[mixing], log = log)
  }
  res
}

# Density using the new parametrisation
dbgev2 = function(x, q, sb, xi, alpha = 0.5, beta = 0.5, p_a = .1, p_b = .2, s = 5, log = FALSE) {
  fix_lengths(x, q, sb, xi, p_a, p_b, s)
  tmp     = new.to.old(c(q,sb,xi), alpha = alpha, beta = beta)
  mu      = tmp$mu
  sigma   = tmp$sigma 
  a       = qgev(p_a, mu, sigma, xi)
  b       = qgev(p_b, mu, sigma, xi)
  res     = rep(NA, length(x))
  gumbel  = which(x <= a)
  frechet = which(x >= b)
  mixing  = which(a < x & x < b)
  if (any(gumbel)) {
    g           = get_gumbel_par(mu[gumbel], sigma[gumbel], xi[gumbel], p_a[gumbel], p_b[gumbel])
    res[gumbel] = dgev(x[gumbel], g$mu, g$sigma, 0, log = log)
  }
  if (any(frechet)) res[frechet] = dgev(x[frechet], mu[frechet], sigma[frechet], xi[frechet], log = log)
  if (any(mixing)) {
    res[mixing] = dbgev_mixing(x[mixing], mu[mixing], sigma[mixing], xi[mixing],
                               p_a[mixing], p_b[mixing], s[mixing], log = log)
  }
  res
}

# Return levels
return_level_bgev = function(period, mu, sigma, xi, p_a = .1, p_b = .2, s = 5) {
  if (any(period <= 1)) warning("invalid period")
  p = ifelse(period > 1, 1 - 1 / period, NA)
  qbgev(p, mu, sigma, xi, p_a, p_b, s)
}

# Tool function for dgev
dbgev_mixing = function(x, mu, sigma, xi, p_a = .1, p_b = .2, s = 5, log = FALSE) {
  g = get_gumbel_par(mu, sigma, xi, p_a, p_b)
  a = qgev(p_a, mu, sigma, xi)
  b = qgev(p_b, mu, sigma, xi)
  if (any(x <= a | x >= b)) stop("x is outside the domain for mixing")
  p = pbeta((x - a) / (b - a), s, s)
  p_der = dbeta((x - a) / (b - a), s, s) / (b - a)
  term1 = - p_der * (1 + xi * (x - mu) / sigma) ^ (-1 / xi)
  term2 = p / sigma * (1 + xi * (x - mu) / sigma) ^ (-1 / xi - 1)
  term3 = p_der * exp(- (x - g$mu) / g$sigma)
  term4 = (1 - p) / g$sigma * exp(- (x - g$mu) / g$sigma)
  term0 = p * log(pgev(x, mu, sigma, xi)) + (1 - p) * log(pgev(x, g$mu, g$sigma, 0))
  res = term0 + log(term1 + term2 + term3 + term4)
  if (!log) res = exp(res)
  res
}

# Tool function for qgev
qbgev_mixing = function(p, mu, sigma, xi, p_a = .1, p_b = .2, s = 5, lower = 0, upper = 100) {
  if (any(p <= p_a | p >= p_b)) stop("p is outside the domain for mixing")
  res = vector("numeric", length(p))
  for (i in seq_along(p)) {
    f = function(x) (pbgev(x, mu, sigma, xi, p_a, p_b, s) - p)[i]
    sol = uniroot(f, lower = lower, upper = upper, extendInt = "upX")
    res[i] = sol$root
  }
  res
}

# Tool function to get parameters for G (Gumbel)
get_gumbel_par = function(mu, sigma, xi, p_a = .1, p_b = .2) {
  if (any(xi < 0)) stop("xi must be nonnegative")
  a = qgev(p_a, mu, sigma, xi)
  b = qgev(p_b, mu, sigma, xi)
  sigma2 = (b - a) / log(log(p_a) / log(p_b))
  mu2 = a + sigma2 * log(-log(p_a))
  list(mu = mu2, sigma = sigma2)
}