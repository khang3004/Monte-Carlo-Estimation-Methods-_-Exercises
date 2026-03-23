#' Exercise 1 - Importance Sampling for E(mu | X = x)
#'
#' @description 
#' Estimates the posterior expectation E(mu | X = x) using Importance Sampling.
#' The data is generated from a Student-t distribution and the proposal is also Student-t.
#'
#' @param n Integer. Number of observations.
#' @param df Numeric. Degrees of freedom for Student-t distribution.
#' @param mu_true Numeric. True location parameter used to generate data.
#' @param scale Numeric. Scale parameter.
#' @param k Integer. Index used for proposal shift (conditioning trick).
#' @param m Integer. Number of importance samples.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return Numeric. Estimated value of E(mu | X = x).
#'

exercise_1 <- function(n = 30, df = 3, mu_true = 0, scale = 1,
                       k = 2, m = 10000, seed = 1) {
  set.seed(seed)
  
  # Generate data
  x <- rt(n, df = df) * scale + mu_true
  
  # Sample mu_j from proposal
  mu_samples <- rt(m, df = df) * scale + x[k]
  
  # Compute weights
  weights <- rep(1, m)
  for (i in 1:n) {
    if (i != k) {
      weights <- weights * (3 + (x[i] - mu_samples)^2)^(-2)
    }
  }
  
  # Estimate expectation
  E_mu_hat <- sum(mu_samples * weights) / sum(weights)
  
  return(E_mu_hat)
}

result <- exercise_1()
cat("Estimated E(mu | X = x):", result, "\n")