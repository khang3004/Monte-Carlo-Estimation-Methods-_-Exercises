#' Importance Sampling for Tail Integral
#' @description Estimates the integral of x^2 * phi(x) from 1 to infinity.
#' @param n_samples Number of simulations (M).
#' @return Estimated integral value and standard error.
solve_importance_sampling <- function(n_samples = 100000) {
  # 1. Generate samples from q(x) = Exp(rate=1) shifted by +1
  # Using rexp(n, rate=1) + 1
  x <- rexp(n_samples, rate = 1) + 1
  
  # 2. Define p(x) - Normal PDF
  p_x <- dnorm(x)
  
  # 3. Define q(x) - Shifted Exponential PDF
  # Note: dexp(x-1, rate=1)
  q_x <- dexp(x - 1, rate = 1)
  
  # 4. Target function f(x) = x^2
  f_x <- x^2
  
  # 5. Calculate Importance Weights and Values
  weights <- p_x / q_x
  estimates <- f_x * weights
  
  # Results
  final_estimate <- mean(estimates)
  std_error <- sd(estimates) / sqrt(n_samples)
  
  return(list(
    estimate = final_estimate,
    se = std_error,
    efficiency = (1/n_samples) / (var(estimates)/n_samples) # Relative efficiency
  ))
}

# Execute
result <- solve_importance_sampling()
print(result[1])
