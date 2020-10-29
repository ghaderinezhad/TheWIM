library(transport)

#' poisson.GammaGamma.vallender
#' 
#' estimated Wasserstein distance between posteriors 
#' from Gamma(a1,b1) vs. Gamma(a2,b2) priors
#' 
#' @param n number of observations
#' @param X sum of observations
#' @param a1 alpha parameter of the reference gamma-prior (shape)
#' @param b1 beta parameter of the reference gamma-prior (rate)
#' @param a2 alpha parameter of the second gamma-prior (shape)
#' @param b2 beta parameter of the second gamma-prior (rate)
#' @return estimated Wasserstein distance 
#' 
poisson.GammaGamma.vallender <- function(n, X, a1, b1, a2, b2) {
  
  # resulting posteriors:
  # p(mu;x) = gamma(shape = X + a, rate = n + b)
  x1 <- rgamma(1000, X + a1, n + b1)  # larger n to get more accurate estimate
  x2 <- rgamma(1000, X + a2, n + b2)
  dw <- wasserstein1d(x1, x2)
  
  return( list(DW = dw, x1 = x1, x2 = x2) )
}
