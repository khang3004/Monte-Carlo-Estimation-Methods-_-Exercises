#' Exercise 3a - Importance Sampling
#' @param N Number of samples
#' @param seed Random seed
#' @return Estimate of E[X^2]

exercise_3a <- function(N = 10000, seed = 123) {
  
  set.seed(seed)
  
  h <- function(x) exp(-abs(x)^3 / 3)
  g <- function(x) dnorm(x)
  
  x <- rnorm(N)
  
  w <- h(x) / g(x)
  EX2_IS <- sum(x^2 * w) / sum(w)
  
  return(EX2_IS)
}

#' Exercise 3b - Riemann Approximation
#' @param N Number of samples
#' @param seed Random seed
#' @return Estimate of E[X^2]

exercise_3b <- function(N = 10000, seed = 123) {
  
  set.seed(seed)
  
  h <- function(x) exp(-abs(x)^3 / 3)
  
  x <- rnorm(N)
  x_sorted <- sort(x)
  
  numerator <- sum(x_sorted[-N]^2 * diff(x_sorted) * h(x_sorted[-N]))
  denominator <- sum(diff(x_sorted) * h(x_sorted[-N]))
  
  EX2_Riemann <- numerator / denominator
  
  return(EX2_Riemann)
}


#' Exercise 3c - Compare IS vs Riemann
#' @param N Number of samples per run
#' @param n_sim Number of simulations
#' @param seed Random seed
#' @return List of summary statistics

exercise_3c <- function(N = 10000, n_sim = 1000, seed = 123) {
  
  set.seed(seed)
  
  results_IS <- numeric(n_sim)
  results_Riemann <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    results_IS[i] <- exercise_3a(N, seed = seed + i)
    results_Riemann[i] <- exercise_3b(N, seed = seed + i)
  }
  
  # Summary statistics
  IS_mean <- mean(results_IS)
  IS_sd   <- sd(results_IS)
  
  R_mean  <- mean(results_Riemann)
  R_sd    <- sd(results_Riemann)
  
  # Print
  cat("Importance Sampling:\n")
  cat("Mean:", IS_mean, " | SD:", IS_sd, "\n\n")
  
  cat("Riemann sum:\n")
  cat("Mean:", R_mean, " | SD:", R_sd, "\n\n")
  
  # Plot
  x_range <- range(c(results_IS, results_Riemann))
  
  hist(results_IS, breaks = 30,
       col = rgb(1, 0, 0, 0.5),
       main = "Comparison: IS vs Riemann",
       xlab = "Estimate",
       xlim = x_range)
  
  hist(results_Riemann, breaks = 30,
       col = rgb(0, 0, 1, 0.5),
       add = TRUE)
  
  legend("topright",
         legend = c("IS", "Riemann"),
         fill = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))
  
  return(list(
    IS_mean = IS_mean,
    IS_sd   = IS_sd,
    Riemann_mean = R_mean,
    Riemann_sd   = R_sd
  ))
}

# a
exercise_3a()

# b
exercise_3b()

# c
result_c <- exercise_3c()
print(result_c)
#Importance Sampling:
  #Mean: 0.7760466  | SD: 0.008364374 

#Riemann sum:
  #Mean: 0.776472  | SD: 5.070328e-05 