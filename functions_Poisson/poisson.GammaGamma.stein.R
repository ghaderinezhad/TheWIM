#' poisson.GammaGamma.stein
#' 
#' bounds on Wasserstein distance between posteriors 
#' from Beta(a,b) (p2) vs. Uniform (p1) priors by Stein's method
#' 
#' @param n number of observations
#' @param X sum of observations
#' @param a1 alpha parameter of the reference gamma-prior (shape)
#' @param b1 beta parameter of the reference gamma-prior (rate)
#' @param a2 alpha parameter of the second gamma-prior (shape)
#' @param b2 beta parameter of the second gamma-prior (rate)
#' @return exact Wasserstein distance (only for a1 < a2 & b1 > b2 or a1 > a2 & b1 < b2)
#' 
poisson.GammaGamma.stein <- function(n, X, a1, b1, a2, b2) {
  
  # Stein's method
  d <- (1/(n + b1))*abs(a2 - a1 - (b2 - b1)*((X + a2)/(n + b2)) ) 
  
  return( d )
} 
