library(transport)

#' poisson.GammaGamma.mopess
#' 
#' calculation of Mean Observed Prior Effective Sample Size (MOPESS) 
#' 
#' @param x original iid dataset
#' @param L maximum feasible MOPESS value
#' @param a1 alpha parameter of the reference gamma prior (shape)
#' @param b1 beta parameter of the reference gamma prior (rate)
#' @param a2 alpha parameter of the second gamma prior (shape)
#' @param b2 beta parameter of the second gamma prior (rate)
#' @return MOPESS
#' 
poisson.GammaGamma.mopess <- function(x, L, a1, b1, a2, b2) {
  
  # replicate process S times
  S     <- 500
  opess <- vector(length = S)
  for (k in 1:S) {
    
    # part a: sampling from posterior predictive distributions
    xl1 <- xl2 <- list()
    mu  <- rgamma(1, sum(x) + a2, length(x) + b2)
    for (m in seq(1, L, 1)) {
      xl1[[m]] <- c(x, rpois(m, mu))
      xl2[[m]] <- c(x, rpois(m, mu))
    }
    
    # part b + c
    W1 <- sapply(xl1, function(xl1) {
      x1 <- rgamma(10000, sum(xl1) + a1, length(xl1) + b1)
      x2 <- rgamma(10000, sum(x) + a2, length(x) + b2)
      dw <- wasserstein1d(x1, x2, p = 2)
      return( dw )
    })
    W2 <- sapply(xl2, function(xl2) {
      x1 <- rgamma(10000, sum(x) + a1, length(x) + b1)
      x2 <- rgamma(10000, sum(xl2) + a2, length(xl2) + b2)
      dw <- wasserstein1d(x1, x2, p = 2)
      return( dw )
    })
    
    # sign function + opess
    s        <- ifelse(min(W1) <= min(W2), 1, -1)
    opess[k] <- ifelse(s == 1, 
                       which.min(W1),
                       -1*which.min(W2))
  }

  return( list(MOPESS = mean(opess), OPESSdist = opess) )
}
  