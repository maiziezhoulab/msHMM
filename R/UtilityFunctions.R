# This page Author: Daniel Lai
# https://bioconductor.org/packages/release/bioc/html/HMMcopy.html

# Normalize a given array to sum to 1
normalize <- function(A) {
  M <- A / (sum(A) + (sum(A) == 0))
  return(M);
}

# Dirichlet probability density function, returns the probability of vector
# x under the Dirichlet distribution with parameter vector alpha
# Author: David Ross http://www.cs.toronto.edu/~dross/code/dirichletpdf.m
dirichletpdf <- function(x, alpha) {
  if (any(x < 0)) {
    return(0);
  }
  if (abs(sum(x) - 1) > 1e-3) {
    stop("Dirichlet PDF: probabilities do not sum to 1")
    return(0);
  }
  p <- exp(lgamma(sum(alpha)) - sum(lgamma(alpha))) * prod(x ^ (alpha - 1))
  return(p)
}

tdistPDF <- function(x, mu, lambda, nu) {
  p <- (gamma(nu / 2 + 0.5)/gamma(nu / 2)) * ((lambda / (pi * nu)) ^ (0.5)) *
    (1 + (lambda * (x - mu) ^ 2) / nu) ^ (-0.5 * nu - 0.5)
  p[is.na(p)] <- 1
  return(p)
}
