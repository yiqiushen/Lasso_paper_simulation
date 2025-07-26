###############################################################################
# 0. Libraries
###############################################################################
library(MASS)    # for mvrnorm()
library(CVXR)    # for Lasso solutions
library(ggplot2)
library(dplyr)
library(gridExtra)

###############################################################################
# 1. Generate Rademacher random vectors
###############################################################################
rsign <- function(n) {
  sample(c(-1, 1), n, replace = TRUE)
}

###############################################################################
# 2. Generate 5×5 block-diagonal X
###############################################################################
generate_block_diagonal_X <- function(n, p, block_size = 5) {
  if (p %% block_size != 0) {
    stop("p must be a multiple of block_size = 5 for this code.")
  }
  
  # Optional coherence
  mu <- 0  # set to 0 for simplicity
  
  # Single 5x5 block: off-diagonal=mu, diag=1
  Sigma_block <- matrix(mu, nrow = block_size, ncol = block_size)
  diag(Sigma_block) <- 1
  
  # Number of blocks
  n_blocks <- p / block_size
  
  # Construct full covariance
  Sigma <- matrix(0, nrow = p, ncol = p)
  for (b in seq_len(n_blocks)) {
    start_col <- (b - 1) * block_size + 1
    end_col   <- b * block_size
    Sigma[start_col:end_col, start_col:end_col] <- Sigma_block
  }
  
  # Generate X ~ N(0, Sigma)
  X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  
  # Scale columns so each has norm sqrt(n)
  col_norms <- sqrt(colSums(X^2))
  scales    <- sqrt(n) / col_norms
  scales[!is.finite(scales)] <- 1
  X <- sweep(X, 2, scales, "*")
  
  return(X)
}

###############################################################################
# 3. Generate data (X, Y, beta, theta)
#    We add outliers in Y by adding sqrt(n)*theta
###############################################################################
generate_data <- function(n, p, s, num_outliers, noise_sd) {
  # (a) Design matrix
  X <- generate_block_diagonal_X(n, p, block_size = 5)
  
  # (b) True beta with s active coordinates
  beta <- rep(0, p)
  beta_support <- sample.int(p, s)
  beta[beta_support] <- runif(s, 8, 16) * rsign(s)
  
  # (c) Outlier vector theta
  theta <- rep(0, n)
  if (num_outliers > 0) {
    outlier_inds <- sample.int(n, num_outliers)
    theta[outlier_inds] <- runif(length(outlier_inds), 16, 32) * rsign(length(outlier_inds))
  }
  
  # (d) Noise
  xi <- rnorm(n, sd = noise_sd)
  
  # (e) Response: Y = X beta + sqrt(n)*theta + xi
  Y <- X %*% beta + sqrt(n)*theta + xi
  
  return(list(X = X, Y = Y, beta = beta, theta = theta))
}

###############################################################################
# 4. Solve Regular Lasso with CVXR:
#    minimize (0.5)*||X * beta - Y||^2 + lambda * ||beta||_1
###############################################################################
solve_lasso_cvxr <- function(X, Y, lambda) {
  n <- nrow(X)
  p <- ncol(X)
  
  beta_var <- Variable(p)
  objective <- 0.5 * sum_squares(X %*% beta_var - Y) + lambda * norm1(beta_var)
  prob <- Problem(Minimize(objective))
  
  result <- solve(prob, solver = "ECOS")  # or "SCS", "OSQP", etc.
  beta_hat <- result$getValue(beta_var)
  return(as.vector(beta_hat))
}

###############################################################################
# 5. Solve Robust Lasso with CVXR, penalizing both beta and outliers:
#    w = (beta, theta) in R^(p+n)
#    X_aug = [X, sqrt(n)*I_n] in R^(n x (p+n))
#    minimize (0.5)*||X_aug * w - Y||^2 + lambda * ||w||_1
#
#    => this penalizes both beta and theta
###############################################################################
solve_augmented_lasso_cvxr <- function(X, Y, n, lambda) {
  p <- ncol(X)
  
  # Build augmented design: [X, sqrt(n)*I_n]
  X_aug <- cbind(X, sqrt(n)*diag(n))
  
  # Single variable w in R^(p + n)
  w <- Variable(p + n)
  
  # L1 penalty on all coefficients => beta + outliers both penalized
  objective <- 0.5 * sum_squares(X_aug %*% w - Y) + lambda * norm1(w)
  prob <- Problem(Minimize(objective))
  
  result <- solve(prob, solver = "ECOS")
  w_hat <- result$getValue(w)
  
  # First p entries of w => beta, last n => outliers
  beta_hat <- w_hat[1:p]
  theta_hat <- w_hat[(p+1):(p+n)]
  
  return(list(beta_hat = as.vector(beta_hat),
              theta_hat = as.vector(theta_hat)))
}

###############################################################################
# 6. Fit both (Regular & Robust) Lasso, compute prediction risk
###############################################################################
fit_and_evaluate <- function(X, Y, beta_true, n, p) {
  # Common lambda
  lambda <- 4 * sqrt(2 * (log(p) + log(n)) / n)
  
  # (A) Solve Regular Lasso
  beta_hat_regular <- solve_lasso_cvxr(X, Y, lambda)
  
  # (B) Solve Robust Lasso (penalize outliers)
  sol_robust <- solve_augmented_lasso_cvxr(X, Y, n, lambda)
  beta_hat_robust <- sol_robust$beta_hat
  # (We won't use sol_robust$theta_hat for the standard "prediction risk" on X*beta)
  
  # Prediction risk = mean((X*(beta_hat - beta_true))^2)
  risk_regular <- mean((X %*% (beta_hat_regular - beta_true))^2)
  risk_robust  <- mean((X %*% (beta_hat_robust  - beta_true))^2)
  
  return(c(regular = risk_regular, robust = risk_robust))
}

###############################################################################
# 7. Main simulation: vary noise_level and outlier_fraction, fix s=15
###############################################################################
set.seed(123)

n <- 1000
p <- 1100
s <- 25
noise_levels   <- c(0.1, 1,2, 5)
outlier_fracs  <- seq(0, 0.15, by = 0.03)
num_repetitions <- 10  # reduce for speed with CVXR

# Create results data frame
results <- expand.grid(
  noise_level      = noise_levels,
  outlier_fraction = outlier_fracs,
  repetition       = 1:num_repetitions
)

results$risk_regular <- NA
results$risk_robust  <- NA

# Loop over all combos
for (i in seq_len(nrow(results))) {
  noise_sd <- results$noise_level[i]
  ofrac    <- results$outlier_fraction[i]
  num_out  <- round(ofrac * n)
  print(c(noise_sd, num_out))
  # Generate data
  data_list <- generate_data(
    n            = n,
    p            = p,
    s            = s,
    num_outliers = num_out,
    noise_sd     = noise_sd
  )
  
  # Fit and evaluate
  risk_vals <- fit_and_evaluate(data_list$X, data_list$Y, data_list$beta, n, p)
  
  results$risk_regular[i] <- risk_vals["regular"]
  results$risk_robust[i]  <- risk_vals["robust"]
}

###############################################################################
# 8. Summarize & plot: 4 lines (Robust) colored by noise level
###############################################################################
summary_results <- results %>%
  group_by(noise_level, outlier_fraction) %>%
  summarise(
    risk_regular = mean(risk_regular),
    risk_robust  = mean(risk_robust),
    .groups      = 'drop'
  )

# (A) Plot robust risk vs outlier fraction, color by noise_level
p_robust <- ggplot(summary_results, aes(x = outlier_fraction,
                                        y = risk_robust,
                                        color = factor(noise_level))) +
  geom_line(size = 1.2) +
  labs(x = "Outlier Fraction",
       y = "Prediction Risk",
       color = "Noise SD") +
  ggtitle(paste0("Augmented Lasso, s=", s, ", n=", n, ", p=", p)) +
  theme_minimal()

# (B) Plot regular risk for comparison
p_regular <- ggplot(summary_results, aes(x = outlier_fraction,
                                         y = risk_regular,
                                         color = factor(noise_level))) +
  geom_line(size = 1.2) +
  labs(x = "Outlier Fraction",
       y = "Prediction Risk",
       color = "Noise SD") +
  ggtitle(paste0("Regular Lasso, s=", s, ", n=", n, ", p=", p)) +
  theme_minimal()

# Display side by side
grid.arrange(p_robust, p_regular, ncol = 2)
