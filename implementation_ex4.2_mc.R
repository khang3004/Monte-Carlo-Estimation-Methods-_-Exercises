#' Estimate Asian Call Option Price using Monte Carlo Simulation
#'
#' @description
#' Estimates the fair value of an Asian Call Option 
#' by simulating daily price paths and calculating the discounted expected payoff.
#'
#' @param S0 Numeric. Current price of the underlying asset at t=0.
#' @param K Numeric. Strike price.
#' @param sigma Numeric. Annual volatility.
#' @param r Numeric. Annual risk-free rate.
#' @param T_days Integer. Time to expiration in days.
#' @param n Integer. Number of simulation paths (default is 100000).
#'
#' @return A numeric value representing the estimated Asian Call Option price.
#' @export
aco_mc_simulation <- function(S0, K, sigma, r, T_days, n = 1e5) {
  # Daily parameters
  t <- T_days / 365
  dt <- 1 / 365
  drift <- (r - 0.5 * sigma^2) * dt
  vol_sqrt_dt <- sigma * sqrt(dt)
  
  # Matrix to store price paths: rows = simulations, cols = time steps (1 to T)
  # S(t+1) = S(t) * exp(drift + vol_sqrt_dt * Z)
  price_paths <- matrix(0, nrow = n, ncol = T_days)
  
  current_S <- rep(S0, n)
  
  for (t in 1:T_days) {
    Z <- rnorm(n)
    current_S <- current_S * exp(drift + vol_sqrt_dt * Z)
    price_paths[, t] <- current_S
  }
  
  # Calculate arithmetic average for each path
  arithmetic_average <- rowMeans(price_paths)
  
  # Calculate discounted payoff: e^(-rT) * max(0, Avg - K)
  payoffs <- pmax(0, arithmetic_average - K)
  discounted_payoffs <- exp(-r * t) * payoffs
  
  return(mean(discounted_payoffs))
}

# Example 4.1 implementation
params <- list(
  S0 = 50, 
  K = 52, 
  sigma = 0.5, 
  r = 0.05, 
  T_days = 30
)

set.seed(123)
asian_price <- aco_mc_simulation(
  params$S0, 
  params$K, 
  params$sigma, 
  params$r, 
  params$T_days
)

cat("Monte Carlo E(A) estimate:", round(asian_price, 4), "\n")