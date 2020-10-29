#' poisson.GammaGamma.neutrality
#' 
#' neutrality of posteriors from Gamma1 (p1) vs. Gamma2 (p2) priors
#' 
#' @param n number of observations
#' @param X sum of observations
#' @param a1 alpha parameter of the reference gamma-prior (shape)
#' @param b1 beta parameter of the reference gamma-prior (rate)
#' @param a2 alpha parameter of the second gamma-prior (shape)
#' @param b2 beta parameter of the second gamma-prior (rate)
#' @return list of neutrality for the two Gamma posteriors
#' 
poisson.GammaGamma.neutrality <- function(n, X, a1, b1, a2, b2) {
  
  # resulting posteriors:
  # p(mu;x) = Gamma(shape = X + a, rate = n + b)
  mle <- X/n
  i1  <- pgamma(mle, shape = X + a1, rate = n + b1)
  i2  <- pgamma(mle, shape = X + a2, rate = n + b2)
  
  return( list(I1 = i1, I2 = i2) )
}


