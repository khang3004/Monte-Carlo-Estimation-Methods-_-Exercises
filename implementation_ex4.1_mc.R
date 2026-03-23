#' Calculate European Call Option Price using Black-Scholes Formula
#'
#' @description
#' Provides the theoretical value for a European Call Option 
#' based on the Black-Scholes model.
#'
#' @param S0 Numeric. Current price of the underlying asset at t=0.
#' @param K Numeric. Strike price of the option.
#' @param sigma Numeric. Annual volatility of the underlying asset's returns.
#' @param r Numeric. Annual risk-free interest rate (e.g., 0.05 for 5%).
#' @param T_days Numeric. Time to expiration expressed in days.
#'
#' @return A numeric value representing the theoretical fair price of the option.
#' @export
eco_black_scholes_formula <- function(S0, K, sigma, r, T_days) {
  t <- T_days / 365
  
  # Calculate d1 and d2 components.
  d1 <- (log(S0 / K) + (r + 0.5 * sigma^2) * t) / (sigma * sqrt(t))
  d2 <- d1 - sigma * sqrt(t)
  
  # Calculate theoretical price using Cumulative Normal Distribution.
  price <- S0 * pnorm(d1) - K * exp(-r * t) * pnorm(d2)
  return(price)
}

#' Estimate European Call Option Price using Monte Carlo Simulation
#'
#' @description
#' Estimates the fair value of a European Call Option by simulating asset price 
#' paths at maturity and calculating the expected discounted payoff.
#'
#' @param S0 Numeric. Current price of the underlying asset at t=0.
#' @param K Numeric. Strike price of the option.
#' @param sigma Numeric. Annual volatility.
#' @param r Numeric. Annual risk-free rate.
#' @param T_days Numeric. Time to expiration in days.
#' @param n Integer. Number of simulation trials (default is 100000).
#'
#' @return A numeric value representing the estimated Monte Carlo price.
#' @export
eco_monte_carlo_simulation <- function(S0, K, sigma, r, T_days, n = 1e5) {
  t <- T_days / 365
  
  # Generate independent standard normal random variables Z ~ N(0, 1).
  Z <- rnorm(n)
  
  # Simulate terminal stock prices S(T) using the formula.
  ST <- S0 * exp((r - 0.5 * sigma^2) * t + sigma * sqrt(t) * Z)
  
  # Calculate discounted payoffs: e^(-rT) * max(0, S(T) - K).
  payoffs <- pmax(0, ST - K)
  discounted_payoffs <- exp(-r * t) * payoffs
  
  # Return the sample mean as the fair price estimate
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

theory_value <- eco_black_scholes_formula(
  params$S0, 
  params$K, 
  params$sigma, 
  params$r, 
  params$T_days)

set.seed(42)
mc_simulation <- eco_monte_carlo_simulation(
  params$S0, 
  params$K, 
  params$sigma, 
  params$r, 
  params$T_days
)

cat("Theoretical fair price:", round(theory_value, 4), "\n")
cat("Monte Carlo fair price estimate:", round(mc_simulation, 4), "\n")
