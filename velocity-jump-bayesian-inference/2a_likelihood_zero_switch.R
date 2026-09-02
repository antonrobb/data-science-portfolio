# MODULE 2: Forward Algorithm and Log-Likelihood

# GENERATIVE AI ATTRIBUTION
#   Passages developed with generative AI assistance are marked inline with
#   [AI-GENERATED] or [AI-ASSISTED] tags. [AI-GENERATED] marks code
#   substantially produced by AI; [AI-ASSISTED] marks code written by me and
#   subsequently debugged, corrected or optimised with AI assistance.
#   Unmarked code is my own. Full details are in the GenAI statement.

# PURPOSE:
#   Compute the log-likelihood of an observed sequence of position
#   increments given a candidate parameter set theta = (v, lambda, sigma).
#   This is the core computational engine that the MCMC sampler (Module 3)
#   will call at every iteration.

# BACKGROUND:
#   In the ion channel setting of Siekmann et al. (2011), the likelihood of
#   an event sequence E given model Q is computed as a chain of matrix
#   multiplications:

#     P(E|Q) = p0 * I_O * A_tau * I_C * A_tau * ... * e

#   where I_O and I_C are projection matrices (binary emission model) and
#   A_tau = exp(Q * tau) is the transition matrix.

#   In this spatial setting, the binary projection matrices are replaced by
#   Gaussian Emission Probabilities. If the particle is in state s during
#   the interval [t_{i-1}, t_i], the expected displacement is v_s * delta_t,
#   and the observed increment delta_y_i is distributed as:

#     delta_y_i | state s ~ N(v_s * delta_t, 2 * sigma^2)

#   (Variance 2*sigma^2, not sigma^2: the increment differences two
#   independent measurement errors. See Section 3 for the full explanation.)

#   The forward algorithm propagates a probability vector alpha_i through
#   time, updating it at each observation using:

#     alpha_i = (alpha_{i-1} * A_tau) * diag(f_i)

#   where f_i is the vector of emission probabilities at step i:
#     sd = sqrt(2) * sigma

#   The log-likelihood is then:
#     log P(increments | theta) = sum of log(normalisation constants)

#   NOTE ON SCALING:
#   The raw forward probabilities become extremely small for long sequences
#   (products of many small numbers underflow to zero). We therefore
#   normalise alpha at each step and accumulate the log of the normalisation
#   constants. This is the standard approach for numerical stability in HMMs.
#   Reference: Rabiner (1989), Sec. V (scaling procedure for the forward variable).

# INPUTS:
#   increments : numeric vector of length N, observed delta_y values
#   v          : numeric vector of length n, velocities per state
#   Q          : n x n generator matrix
#   sigma      : scalar, measurement noise standard deviation
#   delta_t    : scalar, time between observations
#   pi0        : numeric vector of length n, initial state distribution
#                (default: stationary distribution of Q)

# OUTPUT:
#   Scalar log-likelihood value (or -Inf if numerically degenerate)


# SECTION 1: Matrix Exponential via Eigendecomposition

# We need A_tau = exp(Q * delta_t) at every likelihood evaluation.
# R's base expm is not available without a package, so we implement it
# via eigendecomposition:

#   Q = V * diag(eigenvalues) * V^{-1}
#   exp(Q * t) = V * diag(exp(eigenvalues * t)) * V^{-1}

# This is exact for diagonalisable matrices. CTMCs with distinct eigenvalues
# are always diagonalisable, which holds generically.

# Alternatively could install the 'expm' package and use expm::expm(Q * delta_t).
# We use the eigendecomposition here to keep the code self-contained.

compute_transition_matrix <- function(Q, delta_t) {
  n <- nrow(Q)
  
  # Eigendecomposition: Q = V * D * V^{-1}
  # [AI-GENERATED] Matrix exponential by eigendecomposition, including the
  # handling of negligible imaginary parts, clipping of small negative entries
  # and row renormalisation. See GenAI statement.
  eig    <- eigen(Q)
  V      <- eig$vectors          # columns are right eigenvectors
  D      <- eig$values           # eigenvalues (should be real and <= 0 for CTMC)
  V_inv  <- solve(V)             # inverse of eigenvector matrix
  
  # Matrix exponential: exp(Q * delta_t) = V * diag(exp(D * delta_t)) * V^{-1}
  A <- V %*% diag(exp(D * delta_t)) %*% V_inv
  
  # Enforce non-negativity and row-normalisation for numerical stability
  # (tiny imaginary parts and negative values can appear due to floating point)
  A <- Re(A)               # discard negligible imaginary parts
  A <- pmax(A, 0)          # clip any tiny negative values to 0
  A <- A / rowSums(A)      # renormalise rows to sum to 1
  
  return(A)
}


# SECTION 2: Stationary Distribution of Q

# The stationary distribution pi satisfies: pi * Q = 0, sum(pi) = 1.
# This is used as the default initial state distribution pi0.
# We solve the linear system by replacing the last equation of pi * Q = 0
# with the normalisation constraint sum(pi) = 1.

compute_stationary <- function(Q) {
  n <- nrow(Q)
  
  # Build the system: pi * Q = 0 with sum(pi) = 1
  # Transpose Q so we solve Q^T * pi^T = 0
  A_sys <- t(Q)
  
  # Replace last row with normalisation constraint
  A_sys[n, ] <- 1
  
  b <- c(rep(0, n - 1), 1)
  
  # Solve the linear system
  pi_stat <- as.vector(solve(A_sys, b))
  
  # Numerical safety
  pi_stat <- pmax(pi_stat, 0)
  pi_stat <- pi_stat / sum(pi_stat)
  
  return(pi_stat)
}


# SECTION 3: Gaussian Emission Probabilities

# At each time step i, compute the emission probability for each state s:

#   f_i[s] = P(delta_y_i | state = s) = dnorm(delta_y_i, v_s * delta_t,
#                                             sd = sqrt(2) * sigma)

# IMPORTANT — the variance is 2*sigma^2, NOT sigma^2.
# The observed increment is delta_y_i = y_i - y_{i-1} = delta_x_i + delta_eps_i,
# where delta_eps_i = eps_i - eps_{i-1} is the DIFFERENCE of two independent
# N(0, sigma^2) measurement errors. The variance of a difference of two
# independent equal-variance Gaussians is sigma^2 + sigma^2 = 2*sigma^2, so
#   delta_eps_i ~ N(0, 2*sigma^2).
# This matches the data model in Module 1 (noise is added to positions, then
# differenced) and Ceccarelli et al. (2025, Bull. Math. Biol. 87:57), who write
#   delta_y_j ~ N(delta_x_j, 2*sigma^2).
# Using sd = sigma here (the previous version) understates the noise variance
# by a factor of 2 and forces the inferred sigma to inflate to compensate.

# Returns an n-vector of probabilities (one per state).
# These are the spatial analogue of the projection matrices I_O, I_C
# in Siekmann et al. (2011).

compute_emissions <- function(delta_y, v, delta_t, sigma) {
  # Expected displacement in each state
  means <- v * delta_t
  
  # Standard deviation of the differenced noise: sqrt(2) * sigma
  sd_increment <- sqrt(2) * sigma
  
  # Gaussian density evaluated at the observed increment
  emissions <- dnorm(delta_y, mean = means, sd = sd_increment)
  
  return(emissions)  # length-n vector
}


# SECTION 4: The Forward Algorithm (Log-Likelihood)

# Implements the scaled forward algorithm for a hidden Markov model with
# continuous Gaussian emissions and CTMC transitions.

# ALGORITHM (scaled version for numerical stability):

#   Initialisation (i = 0):
#     alpha_0 = pi0 * f_1          (element-wise multiply)
#     c_0     = sum(alpha_0)       (normalisation constant)
#     alpha_0 = alpha_0 / c_0      (normalised forward vector)

#   Recursion (i = 1, ..., N):
#     alpha_i_bar = alpha_{i-1} * A_tau    (propagate through transition matrix)
#     alpha_i     = alpha_i_bar * f_{i+1} (weight by emission probability)
#     c_i         = sum(alpha_i)           (normalisation constant)
#     alpha_i     = alpha_i / c_i          (normalise)

#   Log-likelihood:
#     log P(y | theta) = sum_i log(c_i)

# The normalisation constants c_i accumulate the total probability mass,
# so summing their logs gives the exact log-likelihood.

forward_algorithm <- function(increments, v, Q, sigma, delta_t, pi0 = NULL) {
  
  n <- nrow(Q)   # number of hidden states
  N <- length(increments)
  
  # Compute transition matrix A_tau = exp(Q * delta_t)
  A_tau <- compute_transition_matrix(Q, delta_t)
  
  # Initial state distribution
  if (is.null(pi0)) {
    pi0 <- compute_stationary(Q)
  }
  
  # Storage for log-likelihood accumulation
  log_likelihood <- 0
  
  # Initialisation: weight pi0 by emission at first increment
  f1    <- compute_emissions(increments[1], v, delta_t, sigma)
  alpha <- pi0 * f1           # element-wise: length-n vector
  
  # Scale and accumulate
  c0 <- sum(alpha)
  
  # Guard against numerical zero (degenerate case)
  if (c0 <= 0 || !is.finite(c0)) return(-Inf)
  
  log_likelihood <- log_likelihood + log(c0)
  alpha          <- alpha / c0    # normalised forward vector
  
  # Recursion over remaining increments
  for (i in 2:N) {
    # Step 1: Propagate through transition matrix
    # alpha is a row vector, A_tau is n x n
    # alpha_bar[j] = sum_s alpha[s] * A_tau[s, j]
    alpha_bar <- as.vector(alpha %*% A_tau)
    
    # Step 2: Weight by emission probability at this step
    f_i   <- compute_emissions(increments[i], v, delta_t, sigma)
    alpha <- alpha_bar * f_i
    
    # Step 3: Normalise and accumulate log-likelihood
    c_i <- sum(alpha)
    
    if (c_i <= 0 || !is.finite(c_i)) return(-Inf)
    
    log_likelihood <- log_likelihood + log(c_i)
    alpha          <- alpha / c_i
  }
  
  return(log_likelihood)
}


# SECTION 5: Wrapper — Log-Likelihood from Parameter Vector

# Convenient wrapper used by the MCMC sampler (Module 3).
# Takes a flat parameter vector theta and the data, returns log-likelihood.

# PARAMETER VECTOR CONVENTION (2-state model):
#   theta = c(v1, v2, log_lambda1, log_lambda2, log_sigma)

# We infer log(lambda) and log(sigma) rather than lambda and sigma directly
# because:
#   1. It enforces positivity automatically (exp of anything is positive)
#   2. It treats slow and fast rates symmetrically on a multiplicative scale
#   3. It improves MCMC mixing (the posterior is more Gaussian in log-space)
# This follows Ceccarelli et al. (2025).

# For the 2-state model, P is fixed as [[0,1],[1,0]] so lambda fully
# determines Q. For n > 2 states, P entries would also be in theta.

log_likelihood_from_theta <- function(theta, increments, delta_t, n_states = 2) {
  
  # Unpack theta
  v          <- theta[1:n_states]
  log_lambda <- theta[(n_states + 1):(2 * n_states)]
  log_sigma  <- theta[2 * n_states + 1]
  
  lambda <- exp(log_lambda)
  sigma  <- exp(log_sigma)
  
  # Basic validity checks
  if (sigma <= 0 || any(lambda <= 0)) return(-Inf)
  
  # Build Q (2-state case: P is fixed)
  # For 2 states, the only off-diagonal probabilities are P_12 = P_21 = 1
  if (n_states == 2) {
    P <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  } else {
    stop("n_states > 2 not yet implemented in this wrapper")
  }
  
  Q <- build_Q_matrix(lambda, P)
  
  # Compute and return log-likelihood
  ll <- forward_algorithm(
    increments = increments,
    v          = v,
    Q          = Q,
    sigma      = sigma,
    delta_t    = delta_t,
    pi0        = NULL       # use stationary distribution
  )
  
  return(ll)
}


# SECTION 6: Validation and Testing

# Before using the likelihood in MCMC, we validate it in three ways:

#   Test 1: Log-likelihood at true parameters should be a finite number
#   Test 2: Perturbing parameters should (usually) decrease the likelihood
#   Test 3: Likelihood surface slice along each parameter axis

# These tests give confidence that the forward algorithm is working correctly
# before we hand it to the MCMC sampler.

source("1_simulation.r")   # loads sim object and helper functions

# Run only when this file is executed directly, not when sourced by Modules 2B,
# 2C, 2D or 3.
# [AI-GENERATED] Guard controlling execution when the module is sourced.
if (sys.nframe() == 0) {
  
  cat("\n=== MODULE 2: LIKELIHOOD VALIDATION ===\n\n")
  
  # True parameter vector (in the theta convention used by MCMC)
  # theta = c(v1, v2, log_lambda1, log_lambda2, log_sigma)
  theta_true <- c(
    v_true[1],               # v1 = 2000
    v_true[2],               # v2 = -1500
    log(lambda_true[1]),     # log(lambda1) = log(1.0) = 0
    log(lambda_true[2]),     # log(lambda2) = log(0.5) = -0.693
    log(sigma_true)          # log(sigma)   = log(50)  = 3.912
  )
  
  cat("True theta vector:\n")
  cat("  v1          =", theta_true[1], "\n")
  cat("  v2          =", theta_true[2], "\n")
  cat("  log(lambda1)=", round(theta_true[3], 3), "\n")
  cat("  log(lambda2)=", round(theta_true[4], 3), "\n")
  cat("  log(sigma)  =", round(theta_true[5], 3), "\n\n")
  
  # Test 1: Log-likelihood at true parameters
  ll_true <- log_likelihood_from_theta(theta_true, sim$increments, delta_t)
  cat("Test 1 — Log-likelihood at true parameters:", round(ll_true, 2), "\n")
  cat("  (Should be a finite negative number)\n\n")
  
  # Test 2: Perturbed parameters should give lower likelihood
  set.seed(1)
  
  # Scale perturbations relative to each parameter's magnitude
  perturb_sd <- c(200, 200, 0.3, 0.3, 0.3)
  
  # Number of random perturbations to test, and storage for their log-likelihoods
  n_perturb  <- 20
  ll_perturb <- numeric(n_perturb)
  
  for (k in 1:n_perturb) {
    theta_perturb <- theta_true + rnorm(length(theta_true), sd = perturb_sd)
    ll_perturb[k] <- log_likelihood_from_theta(theta_perturb, sim$increments, delta_t)
  }
  
  cat("Test 2 — Log-likelihoods at", n_perturb, "randomly perturbed parameter sets:\n")
  cat("  True LL  :", round(ll_true, 2), "\n")
  cat("  Mean pert:", round(mean(ll_perturb), 2), "\n")
  cat("  Max pert :", round(max(ll_perturb), 2), "\n")
  cat("  Fraction below true:", mean(ll_perturb < ll_true), "\n\n")
  
  # Test 3: Likelihood surface slice along each parameter
  cat("Test 3 — Likelihood surface slices (plotting)...\n")
  
  param_names  <- c("v1", "v2", "log(lambda1)", "log(lambda2)", "log(sigma)")
  param_ranges <- list(
    seq(1500,  2500, length.out = 40),   # v1
    seq(-2000, -800, length.out = 40),   # v2
    seq(-1.5,   1.5, length.out = 40),   # log(lambda1)
    seq(-2.0,   0.5, length.out = 40),   # log(lambda2)
    seq(2.5,    6.0, length.out = 40)    # log(sigma)
  )
  
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  
  for (p in 1:5) {
    param_grid <- param_ranges[[p]]
    ll_grid    <- numeric(length(param_grid))
    
    for (k in seq_along(param_grid)) {
      theta_k    <- theta_true
      theta_k[p] <- param_grid[k]
      ll_grid[k] <- log_likelihood_from_theta(theta_k, sim$increments, delta_t)
    }
    
    plot(param_grid, ll_grid,
         type = "l", lwd = 2,
         col  = "dodgerblue3",
         xlab = param_names[p],
         ylab = "Log-likelihood",
         main = paste("Slice along", param_names[p]))
    
    # Mark true value
    abline(v   = theta_true[p], col = "orange2", lty = 2, lwd = 2)
    legend("bottomright",
           legend = c("LL slice", "True value"),
           col    = c("dodgerblue3", "orange2"),
           lty    = c(1, 2), lwd = 2, bty = "n", cex = 0.8)
  }
  
  par(mfrow = c(1, 1))
  
  # 2D joint slice: v1 vs log(sigma)
  v1_grid    <- seq(1500, 2500, length.out = 30)
  sigma_grid <- seq(2.5, 6.0, length.out = 30)
  ll_2d      <- matrix(NA, nrow = length(v1_grid), ncol = length(sigma_grid))
  
  for (i in seq_along(v1_grid)) {
    for (j in seq_along(sigma_grid)) {
      theta_ij    <- theta_true
      theta_ij[1] <- v1_grid[i]
      theta_ij[5] <- sigma_grid[j]
      ll_2d[i, j] <- log_likelihood_from_theta(theta_ij, sim$increments, delta_t)
    }
  }
  
  contour(v1_grid, sigma_grid, ll_2d,
          nlevels = 20,
          xlab    = "v1",
          ylab    = "log(sigma)",
          main    = "Joint likelihood surface: v1 vs log(sigma)")
  
  # Mark true values
  abline(v = 2000, col = "orange2", lty = 2)
  abline(h = log(sigma_true), col = "orange2", lty = 2)
}