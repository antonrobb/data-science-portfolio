# MODULE 3: Bayesian Inference via Adaptive Metropolis-Hastings MCMC

# GENERATIVE AI ATTRIBUTION
#   Passages developed with generative AI assistance are marked inline with
#   [AI-GENERATED] or [AI-ASSISTED] tags. [AI-GENERATED] marks code
#   substantially produced by AI; [AI-ASSISTED] marks code written by me and
#   subsequently debugged, corrected or optimised with AI assistance.
#   Unmarked code is my own. Full details are in the GenAI statement.

# PURPOSE:
#   Estimate the full posterior distribution of the parameter vector
#   theta = (v1, v2, log_lambda1, log_lambda2, log_sigma) given the observed
#   position increments from Module 1, using the exact forward-algorithm
#   log-likelihood from Module 2.

#   Where a maximum-likelihood approach returns a single point estimate, the
#   Bayesian approach returns a full distribution for every parameter,
#   quantifying uncertainty, which is the motivation for using a Bayesian
#   rather than a point-estimate approach to infer the hidden dynamics.

# BACKGROUND:
#   The posterior is, by Bayes' theorem:

#     p(theta | data) propto  L(data | theta) * p(theta)

#   We work entirely on the LOG scale to avoid underflow:

#     log p(theta | data) = log L(data | theta) + log p(theta) + const

#   The log-likelihood log L(data | theta) is provided by Module 2's
#   log_likelihood_from_theta(). This module supplies the log-prior and the
#   Metropolis-Hastings machinery that explores the posterior.

#   ALGORITHM (Adaptive Metropolis-Hastings; Ceccarelli et al.
#   2025 (Bull. Math. Biol. 87:57), Siekmann et al. 2011, Haario et al. 2001):
#     1. Start at theta_0 (drawn from the prior or set near a sensible guess).
#     2. Propose theta* ~ Normal(theta_current, Sigma), a random-walk proposal.
#     3. Compute the log acceptance ratio
#          log r = [logpost(theta*)] - [logpost(theta_current)]
#     4. Accept theta* with probability min(1, exp(log r)); else stay.
#     5. During burn-in, periodically adapt Sigma to the empirical covariance
#        of the chain so far (scaled by 2.38^2 / d, the optimal factor of
#        Roberts & Rosenthal 2001), which tunes the acceptance rate toward
#        the ~0.234 optimum for random-walk Metropolis.
#     6. Run several independent chains from dispersed starts and assess
#        convergence with the Gelman-Rubin statistic R-hat.

#   The random-walk proposal is symmetric, so the proposal density cancels in
#   the acceptance ratio (this is plain Metropolis, not full
#   Metropolis-Hastings with an asymmetric correction term).

# INPUTS (to run_mcmc):
#   increments : numeric vector of observed delta_y (from Module 1)
#   delta_t    : time between observations
#   n_states   : number of hidden states (2 for now)
#   prior      : a list describing the prior (see make_prior below)
#   n_iter     : iterations per chain
#   burn_in    : number of initial iterations to discard / adapt over
#   ... plus adaptation and initialisation controls

# OUTPUTS:
#   A list with the sampled chains, acceptance rates, log-posterior trace,
#   and the parameter names, ready for Module 4 (diagnostics).


source("1_simulation.r")   # simulation + build_Q_matrix
source("2a_likelihood_zero_switch.r")    # Option A: log_likelihood_from_theta
source("2b_likelihood_one_switch.r")   # Option B: log_likelihood_from_theta_B
source("2c_likelihood_two_switch.r")   # Option C: log_likelihood_from_theta_C
source("2d_likelihood_exact.r")   # Option D: log_likelihood_from_theta_D


# SECTION 1: Prior Distribution

# We place independent priors on each component of theta.

#   v1, v2  (velocities, real-valued):
#       Normal(mean = 0, sd = v_prior_sd), a weakly-informative prior. It is
#       centred at zero so it does not bias the sign of the velocity, but its
#       width keeps velocities in a physically plausible range. Ceccarelli et
#       al. use uniform priors; a broad Normal is the smooth analogue and
#       mixes better. Set v_prior_sd large (e.g. 5000) to stay uninformative.

#   log_lambda1, log_lambda2  (log switching rates):
#       Normal(mean = log_lambda_mean, sd = log_lambda_sd) on the LOG scale.
#       A Normal prior on log(lambda) is a log-Normal prior on lambda: it
#       enforces positivity, treats fast and slow rates symmetrically on a
#       multiplicative scale, and gently penalises extreme rates. This is the
#       smooth counterpart of Siekmann's exponential penalty on rate
#       constants.

#   log_sigma  (log measurement-noise sd):
#       Normal(mean = log_sigma_mean, sd = log_sigma_sd) on the LOG scale,
#       i.e. a log-Normal prior on sigma. Again this enforces sigma > 0 and is
#       weakly informative.

# Working on the log scale for lambda and sigma means the prior is evaluated
# directly in the same coordinates theta lives in, so no Jacobian term is
# needed: inference is carried out directly in the transformed space.

make_prior <- function(n_states      = 2,
                       v_prior_sd     = 5000,
                       log_lambda_mean = 0,    # prior median rate = exp(0) = 1
                       log_lambda_sd   = 2,    # broad: covers ~0.02 to ~50
                       log_sigma_mean  = log(50),
                       log_sigma_sd    = 1) {  # broad around sigma = 50
  list(
    n_states        = n_states,
    v_prior_sd      = v_prior_sd,
    log_lambda_mean = log_lambda_mean,
    log_lambda_sd   = log_lambda_sd,
    log_sigma_mean  = log_sigma_mean,
    log_sigma_sd    = log_sigma_sd
  )
}


# Evaluate the log-prior density at a parameter vector theta.
# theta = c(v_1, ..., v_n, log_lambda_1, ..., log_lambda_n, log_sigma)
log_prior <- function(theta, prior) {
  n <- prior$n_states
  
  v          <- theta[1:n]
  log_lambda <- theta[(n + 1):(2 * n)]
  log_sigma  <- theta[2 * n + 1]
  
  # IDENTIFIABILITY CONSTRAINT (label-switching fix).
  # The likelihood is invariant under permuting the state labels: swapping
  # (v_i, lambda_i) with (v_j, lambda_j) gives an identical likelihood, so the
  # posterior has n! symmetric modes. With unconstrained chains, some chains
  # settle on one labelling and others on a permuted one, which inflates the
  # Gelman-Rubin R-hat for the velocities even though every chain has found
  # the same physical solution.
  
  # We break the symmetry by requiring the velocities to be ordered, v_1 >
  # v_2 > ... > v_n. Any proposal violating the order is assigned zero prior
  # density (-Inf log-prior), so all chains are confined to the single
  # canonical labelling. This is the standard remedy for label switching in
  # mixture and hidden-state models (e.g. Stephens 2000; Jasra et al. 2005).
  # [AI-GENERATED] Compact test for the velocity ordering constraint.
  if (is.unsorted(rev(v), strictly = TRUE)) {
    return(-Inf)   # v is not strictly decreasing -> outside the ordered region
  }
  
  lp <- 0
  
  # Velocities: independent Normal(0, v_prior_sd), restricted to v_1>...>v_n.
  # (The truncation to the ordered region multiplies the density by a constant
  # n!, which is irrelevant for MCMC since it cancels in every acceptance
  # ratio; we therefore omit it.)
  lp <- lp + sum(dnorm(v, mean = 0, sd = prior$v_prior_sd, log = TRUE))
  
  # Log switching rates: independent Normal(log_lambda_mean, log_lambda_sd)
  lp <- lp + sum(dnorm(log_lambda,
                       mean = prior$log_lambda_mean,
                       sd   = prior$log_lambda_sd, log = TRUE))
  
  # Log noise sd: Normal(log_sigma_mean, log_sigma_sd)
  lp <- lp + dnorm(log_sigma,
                   mean = prior$log_sigma_mean,
                   sd   = prior$log_sigma_sd, log = TRUE)
  
  return(lp)
}


# SECTION 2: Log-Posterior

# The (unnormalised) log-posterior is log-likelihood + log-prior.
# We add a guard: if either term is non-finite (e.g. degenerate parameters
# that made the forward algorithm return -Inf), the whole log-posterior is
# -Inf so the proposal is rejected.

# The `likelihood` argument selects which model's likelihood to use:
#   "A" -> log_likelihood_from_theta   (up-to-zero-switch; Module 2)
#   "B" -> log_likelihood_from_theta_B (up-to-one-switch;  Module 2B)
#   "C" -> log_likelihood_from_theta_C (up-to-two-switch;  Module 2C)
#   "D" -> log_likelihood_from_theta_D (exact;             Module 2D)
# This lets the same sampler produce all four posteriors for a like-for-like
# comparison; only this one argument changes between runs.

# `fix_sigma`: if given a numeric value, sigma is held at that value and is NOT
# sampled. In that case `theta` has length 2n (velocities and log-rates only)
# and log(fix_sigma) is appended before evaluating the prior and likelihood.
# Note that constraining sigma with a very tight prior is NOT equivalent: in
# testing it consistently destabilised the sampler (R-hat 1.35 to 3.88), so
# sigma must be removed from the parameter vector rather than merely pinned.

log_posterior <- function(theta, increments, delta_t, prior,
                          likelihood = "A", fix_sigma = NULL) {
  if (!is.null(fix_sigma)) theta <- c(theta, log(fix_sigma))
  lp <- log_prior(theta, prior)
  if (!is.finite(lp)) return(-Inf)
  
  ll <- switch(likelihood,
               "A" = log_likelihood_from_theta(theta, increments, delta_t,
                                               n_states = prior$n_states),
               "B" = log_likelihood_from_theta_B(theta, increments, delta_t,
                                                 n_states = prior$n_states),
               "C" = log_likelihood_from_theta_C(theta, increments, delta_t,
                                                 n_states = prior$n_states),
               "D" = log_likelihood_from_theta_D(theta, increments, delta_t,
                                                 n_states = prior$n_states),
               stop("likelihood must be 'A', 'B', 'C', or 'D'")
  )
  if (!is.finite(ll)) return(-Inf)
  
  return(ll + lp)
}


# SECTION 3: Single-Chain Adaptive Metropolis-Hastings

# Runs one chain of length n_iter. During the burn-in phase the proposal
# covariance Sigma is periodically re-estimated from the chain history,
# following Haario et al. (2001) and the optimal scaling of Roberts &
# Rosenthal (2001): Sigma = (2.38^2 / d) * cov(samples) + epsilon * I.
# A small epsilon * I term keeps Sigma positive-definite and prevents the
# proposal from collapsing.

run_chain <- function(theta_init,
                      increments,
                      delta_t,
                      prior,
                      n_iter        = 10000,
                      burn_in       = 2000,
                      Sigma_init    = NULL,
                      adapt_every   = 100,
                      adapt_start   = 300,   # iteration at which adaptation begins
                      # (default: halfway through burn-in)
                      epsilon       = 1e-6,
                      likelihood    = "A",   # "A"/"B"/"C"/"D"
                      fix_sigma     = NULL,  # if numeric, sigma is not sampled
                      v_floor       = 100,   # proposal SD floor for velocities
                      verbose       = TRUE) {
  
  d <- length(theta_init)
  
  # Optimal random-walk scaling factor (Roberts & Rosenthal 2001)
  scale_factor <- (2.38^2) / d
  
  # Begin adapting only once the chain has had time to migrate toward the
  # mode. Adapting on the early transient estimates the covariance of a
  # moving chain, not the posterior, and can lock the sampler in place.
  if (is.null(adapt_start)) adapt_start <- floor(burn_in / 2)
  
  # Floor on the per-parameter proposal SD. Prevents a degenerate (near-zero
  # variance) covariance estimate from collapsing the proposal so that the
  # chain freezes. Matched to the magnitude of each parameter.
  n_st       <- prior$n_states
  # [AI-ASSISTED] Lower bound on the proposal standard deviation. The value
  # required differs by sample size: a bound of 20 was too small at N = 200 and
  # the velocity chains drifted, while the bound of 100 that resolved this was
  # roughly seven times too wide at N = 1600 and the acceptance rate collapsed
  # to between 0.001 and 0.03. Both failures were diagnosed with AI assistance.
  # See GenAI statement and Section 4.11 of the report.
  scale_floor <- c(rep(v_floor^2, n_st), # velocities: SD floor (v_floor)
                   rep(0.01^2, n_st))  # log-lambda: SD floor 0.01
  if (is.null(fix_sigma)) scale_floor <- c(scale_floor, 0.01^2)  # log-sigma
  
  # Initial proposal covariance: if none supplied, use a modest diagonal.
  # The scales reflect the rough magnitude of each parameter:
  #   velocities are on the scale of (increment / delta_t), i.e. hundreds-to-
  #   thousands, so a proposal SD of ~300 lets a chain traverse that range in
  #   a reasonable number of steps; log-rates and log-sigma are order 1.
  if (is.null(Sigma_init)) {
    n <- prior$n_states
    diag_scales <- c(rep(400^2, n),      # velocity variances (SD = 400)
                     rep(0.1^2, n))      # log-lambda variances
    if (is.null(fix_sigma)) diag_scales <- c(diag_scales, 0.1^2)  # log-sigma
    Sigma <- diag(diag_scales, nrow = d)
  } else {
    Sigma <- Sigma_init
  }
  
  # Storage
  chain      <- matrix(NA_real_, nrow = n_iter, ncol = d)
  logpost    <- numeric(n_iter)
  accept_vec <- logical(n_iter)
  
  # Initialise
  theta_curr <- theta_init
  lp_curr    <- log_posterior(theta_curr, increments, delta_t, prior, likelihood, fix_sigma)
  
  # If the starting point is impossible, nudge until it is finite
  tries <- 0
  while (!is.finite(lp_curr) && tries < 100) {
    theta_curr <- theta_init + rnorm(d, sd = c(rep(50, prior$n_states),
                                               rep(0.1, prior$n_states), 0.1))
    lp_curr <- log_posterior(theta_curr, increments, delta_t, prior, likelihood, fix_sigma)
    tries   <- tries + 1
  }
  if (!is.finite(lp_curr)) {
    stop("Could not find a finite-posterior starting point. Check theta_init / prior.")
  }
  
  # Cholesky factor of Sigma for efficient multivariate-normal proposals
  # [AI-GENERATED] begin: Cholesky factorisation of the proposal covariance,
  # computed once outside the sampling loop, so that multivariate normal
  # proposals are drawn as theta + crossprod(chol(Sigma), z). See GenAI statement.
  chol_Sigma <- chol(Sigma)
  
  for (it in 1:n_iter) {
    # Propose: theta* = theta_curr + L^T z,  z ~ N(0, I)
    z         <- rnorm(d)
    theta_star <- theta_curr + as.vector(crossprod(chol_Sigma, z))
    # [AI-GENERATED] end
    
    lp_star <- log_posterior(theta_star, increments, delta_t, prior, likelihood, fix_sigma)
    
    # Metropolis acceptance (symmetric proposal: no q-ratio correction)
    log_r <- lp_star - lp_curr
    if (is.finite(log_r) && log(runif(1)) < log_r) {
      theta_curr <- theta_star
      lp_curr    <- lp_star
      accept_vec[it] <- TRUE
    }
    
    chain[it, ] <- theta_curr
    logpost[it] <- lp_curr
    
    # Adaptation during burn-in
    # Only adapt after adapt_start, and estimate the covariance from the
    # post-adapt_start window (the chain's behaviour once it has settled),
    # not from the full transient history.
    # [AI-GENERATED] begin: adaptive covariance update. Rolling-window covariance,
    # the 2.38^2/d optimal scaling of Roberts and Rosenthal (2001), the ridge term
    # maintaining positive definiteness, the lower bound on the diagonal, and the
    # recovery step retaining the previous covariance if factorisation fails.
    if (it <= burn_in && it >= adapt_start && it %% adapt_every == 0) {
      window_idx <- adapt_start:it
      if (length(window_idx) > d + 1) {
        emp_cov <- cov(chain[window_idx, , drop = FALSE])
        Sigma   <- scale_factor * emp_cov + epsilon * diag(d)
        
        # Enforce the proposal-SD floor on the diagonal so the chain can
        # never freeze on a degenerate covariance estimate.
        dd <- diag(Sigma)
        diag(Sigma) <- pmax(dd, scale_floor)
        
        # Refresh Cholesky; fall back to previous Sigma if not PD
        chol_try <- tryCatch(chol(Sigma), error = function(e) NULL)
        if (!is.null(chol_try)) chol_Sigma <- chol_try
      }
      # [AI-GENERATED] end
      if (verbose) {
        cat(sprintf("  [adapt] iter %5d  acc(last %d) = %.3f\n",
                    it, adapt_every,
                    mean(accept_vec[(it - adapt_every + 1):it])))
      }
    }
  }
  
  acc_rate_overall <- mean(accept_vec)
  acc_rate_post    <- mean(accept_vec[(burn_in + 1):n_iter])
  
  return(list(
    chain        = chain,        # n_iter x d matrix of samples
    logpost      = logpost,      # log-posterior trace
    accept       = accept_vec,   # per-iteration accept flags
    acc_rate     = acc_rate_overall,
    acc_rate_post = acc_rate_post,
    Sigma_final  = Sigma,        # final adapted proposal covariance
    burn_in      = burn_in,
    n_iter       = n_iter
  ))
}


# SECTION 4: Multi-Chain Driver

# Runs n_chains independent chains from dispersed starting points. Dispersed
# starts are essential for the Gelman-Rubin diagnostic to be meaningful: if
# chains started in different parts of parameter space converge to the same
# distribution, that is strong evidence of convergence.

run_mcmc <- function(increments,
                     delta_t,
                     prior        = make_prior(),
                     n_chains     = 4,
                     n_iter       = 10000,
                     burn_in      = 2000,
                     theta_centre = NULL,
                     start_disp   = NULL,
                     seed         = NULL,
                     likelihood    = "A",   # "A"/"B"/"C"/"D"
                     fix_sigma     = NULL,  # if numeric, sigma is not sampled
                     v_floor      = 100,   # proposal SD floor for velocities
                     Sigma_init   = NULL,  # initial proposal covariance (d x d)
                     verbose       = TRUE) {
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- prior$n_states
  d <- if (is.null(fix_sigma)) 2 * n + 1 else 2 * n
  
  # Default centre for dispersing starts. The data fixes the velocity scale
  # for us: increments / delta_t are observed displacements per unit time, so
  # their spread brackets the plausible velocity range. We centre velocities
  # at +/- this scale (one positive, one negative state) rather than at 0, so
  # chains start near the data-implied velocity scale rather than thousands of
  # units away from it. Rates and sigma start at their prior means.
  v_scale <- stats::sd(increments) / delta_t   # rough velocity magnitude
  if (is.null(theta_centre)) {
    v_centre <- v_scale * rep(c(1, -1), length.out = n)  # alternate signs
    theta_centre <- c(v_centre,
                      rep(prior$log_lambda_mean, n))   # log-rates at prior mean
    if (is.null(fix_sigma)) {
      theta_centre <- c(theta_centre, prior$log_sigma_mean)  # log-sigma
    }
  }
  
  # Per-parameter dispersion of the starting points across chains. Velocities
  # are dispersed on the same scale as the data so starts remain overdispersed
  # (good for Gelman-Rubin) without being absurdly far from the mode.
  if (is.null(start_disp)) {  
    start_disp <- c(rep(v_scale, n),  # velocities: spread ~ data velocity scale
                    rep(1.0, n))    # log-lambda: spread on log scale
    if (is.null(fix_sigma)) start_disp <- c(start_disp, 0.5)  # log-sigma
  }
  chains <- vector("list", n_chains)
  
  for (c in 1:n_chains) {
    theta_init <- theta_centre + rnorm(d, sd = start_disp)
    
    # Enforce the velocity ordering v_1 > ... > v_n on the starting point so
    # it lies inside the constrained prior support (otherwise its prior is
    # -Inf and the chain wastes the nudge loop).
    theta_init[1:n] <- sort(theta_init[1:n], decreasing = TRUE)
    
    if (verbose) {
      cat(sprintf("\n=== Chain %d / %d ===\n", c, n_chains))
      cat("  start:", paste(round(theta_init, 2), collapse = ", "), "\n")
    }
    
    chains[[c]] <- run_chain(
      theta_init = theta_init,
      increments = increments,
      delta_t    = delta_t,
      prior      = prior,
      n_iter     = n_iter,
      burn_in    = burn_in,
      likelihood = likelihood,
      fix_sigma  = fix_sigma,
      v_floor    = v_floor,
      Sigma_init = Sigma_init,
      verbose    = verbose
    )
    
    if (verbose) {
      cat(sprintf("  post-burn-in acceptance rate: %.3f\n",
                  chains[[c]]$acc_rate_post))
    }
  }
  
  param_names <- c(paste0("v", 1:n),
                   paste0("log_lambda", 1:n))
  if (is.null(fix_sigma)) param_names <- c(param_names, "log_sigma")
  
  return(list(
    chains      = chains,
    param_names = param_names,
    n_chains    = n_chains,
    n_iter      = n_iter,
    burn_in     = burn_in,
    prior       = prior,
    delta_t     = delta_t,
    likelihood  = likelihood,
    fix_sigma   = fix_sigma
  ))
}


# SECTION 5: Gelman-Rubin Convergence Diagnostic (R-hat)

# For each parameter, R-hat compares the variance between chains to the
# variance within chains. If the chains have converged to the same
# distribution, between-chain and within-chain variance agree and R-hat -> 1.
# If R-hat < 1.1 it indicates acceptable convergence (Gelman & Rubin 1992).
# We use the post-burn-in samples only.

# [AI-GENERATED] Gelman-Rubin statistic implemented directly from its
# definition rather than through an existing package.
gelman_rubin <- function(mcmc_out) {
  chains  <- mcmc_out$chains
  burn_in <- mcmc_out$burn_in
  n_iter  <- mcmc_out$n_iter
  m       <- mcmc_out$n_chains            # number of chains
  d       <- length(mcmc_out$param_names) # number of parameters
  
  # Post-burn-in length per chain
  post_idx <- (burn_in + 1):n_iter
  n        <- length(post_idx)
  
  rhat <- numeric(d)
  
  for (p in 1:d) {
    # Extract the post-burn-in samples for parameter p from each chain:
    # a matrix with n rows (samples) and m columns (chains)
    M <- sapply(chains, function(ch) ch$chain[post_idx, p])
    
    chain_means <- colMeans(M)        # mean of each chain
    grand_mean  <- mean(chain_means)
    
    # Between-chain variance B and within-chain variance W
    B <- (n / (m - 1)) * sum((chain_means - grand_mean)^2)
    W <- mean(apply(M, 2, var))
    
    # Marginal posterior variance estimate (pooled)
    var_hat <- ((n - 1) / n) * W + (1 / n) * B
    
    rhat[p] <- sqrt(var_hat / W)
  }
  
  names(rhat) <- mcmc_out$param_names
  return(rhat)
}


# SECTION 6: Posterior Summary

# Pools the post-burn-in samples from all chains and reports, for each
# parameter, the posterior mean, sd, and 95% credible interval. Switching
# rates and sigma are also reported on their natural (exponentiated) scale.

posterior_summary <- function(mcmc_out) {
  chains  <- mcmc_out$chains
  burn_in <- mcmc_out$burn_in
  n_iter  <- mcmc_out$n_iter
  names_p <- mcmc_out$param_names
  n       <- mcmc_out$prior$n_states
  
  post_idx <- (burn_in + 1):n_iter
  
  # Pool all chains: stack post-burn-in samples into one big matrix
  pooled <- do.call(rbind, lapply(chains, function(ch) ch$chain[post_idx, , drop = FALSE]))
  colnames(pooled) <- names_p
  
  # Summary on the sampling (theta) scale
  summ <- data.frame(
    parameter = names_p,
    mean      = apply(pooled, 2, mean),
    sd        = apply(pooled, 2, sd),
    q2.5      = apply(pooled, 2, quantile, probs = 0.025),
    q50       = apply(pooled, 2, quantile, probs = 0.5),
    q97.5     = apply(pooled, 2, quantile, probs = 0.975),
    row.names = NULL
  )
  
  # Add natural-scale summaries for lambda and sigma. When sigma was held
  # fixed it is not part of the parameter vector, so only the rates are
  # exponentiated.
  lambda_cols <- (n + 1):(2 * n)
  exp_cols    <- lambda_cols
  if (is.null(mcmc_out$fix_sigma)) exp_cols <- c(exp_cols, 2 * n + 1)
  
  nat <- pooled
  nat[, exp_cols] <- exp(pooled[, exp_cols, drop = FALSE])
  nat_names <- names_p
  nat_names[lambda_cols] <- paste0("lambda", 1:n)
  if (is.null(mcmc_out$fix_sigma)) nat_names[2 * n + 1] <- "sigma"
  
  summ_nat <- data.frame(
    parameter = nat_names,
    mean      = apply(nat, 2, mean),
    sd        = apply(nat, 2, sd),
    q2.5      = apply(nat, 2, quantile, probs = 0.025),
    q50       = apply(nat, 2, quantile, probs = 0.5),
    q97.5     = apply(nat, 2, quantile, probs = 0.975),
    row.names = NULL
  )
  
  return(list(
    theta_scale   = summ,       # on (v, log_lambda, log_sigma) scale
    natural_scale = summ_nat,   # on (v, lambda, sigma) scale
    pooled        = pooled      # pooled samples, for plotting in Module 4
  ))
}


# SECTION 7: Run the Sampler — Options A, B, C and D compared

# Runs the same sampler, same data, same prior, same seed under both
# likelihoods, so the only difference between the two posteriors is the
# within-interval switching treatment. Option B should recover sigma much 
# closer to the truth (~50) and wider, more accurate velocities than the
# biased Option A (sigma ~105).

cat("\n=== MODULE 3: MCMC INFERENCE (Options A, B, C, D) ===\n")

prior <- make_prior(
  n_states        = 2,
  v_prior_sd      = 5000,        # weakly-informative on velocities
  log_lambda_mean = 0,
  log_lambda_sd   = 2,
  log_sigma_mean  = log(50),
  log_sigma_sd    = 1
)

# Sampler settings are shared across options so the comparison is like-for-like,
# but they are exposed as arguments because Option C needs a longer run: on the
# shared settings its log_lambda2 R-hat was 1.115, above the 1.1 threshold, so
# it is run with six chains of 12000 iterations (4000 burn-in) to bring all
# parameters below 1.01. The change is a response to a convergence diagnostic,
# not a change of model, and the longer run leaves Option C with more retained
# samples than the others; this affects only the Monte-Carlo precision of its
# quantiles, which is already negligible at 16000 samples.
run_one <- function(which_lik, n_chains = 4, n_iter = 6000, burn_in = 2000) {
  cat(sprintf("\n########## OPTION %s ##########\n", which_lik))
  out <- run_mcmc(
    increments = sim$increments,
    delta_t    = delta_t,
    prior      = prior,
    n_chains   = n_chains,
    n_iter     = n_iter,
    burn_in    = burn_in,
    seed       = 123,           # same seed -> identical starts for fair compare
    likelihood = which_lik,
    verbose    = TRUE
  )
  cat(sprintf("\n--- Option %s: Gelman-Rubin R-hat (want < 1.1) ---\n", which_lik))
  print(round(gelman_rubin(out), 4))
  cat(sprintf("\n--- Option %s: Posterior summary (natural scale) ---\n", which_lik))
  print(posterior_summary(out)$natural_scale, digits = 4)
  out
}

mcmc_A <- run_one("A")
mcmc_B <- run_one("B")
mcmc_C <- run_one("C", n_chains = 6, n_iter = 12000, burn_in = 4000)
mcmc_D <- run_one("D")

# Keep mcmc_out pointing at Option B by default for Module 4
mcmc_out <- mcmc_B

cat("\n--- True values for comparison ---\n")
cat("  v1 =", v_true[1], "  v2 =", v_true[2], "\n")
cat("  lambda1 =", lambda_true[1], "  lambda2 =", lambda_true[2], "\n")
cat("  sigma =", sigma_true, "\n")

saveRDS(mcmc_A, "mcmc_A.rds")
saveRDS(mcmc_B, "mcmc_B.rds")
saveRDS(mcmc_C, "mcmc_C.rds")
saveRDS(mcmc_D, "mcmc_D.rds")

# SECTION 8: Experiment — faster switching

# Options A/B/C/D are near-indistinguishable at the baseline rates because
# only about 0.3 switches occur per observation interval, so truncating the
# switch count costs almost nothing. This experiment regenerates the test data
# with faster switching, holding v, sigma, delta_t and T fixed, so that
# intervals containing two or more switches become common and the truncated
# likelihoods are put under strain. lambda2 = lambda1 / 2 throughout, matching
# the sweep in Ceccarelli et al. (2025).

# The log-lambda prior is widened to sd = 3 for this experiment: the default
# sd = 2 is centred on lambda = 1 and would fight a true rate of 4.

run_rate_experiment <- function(lambda1, likelihoods = c("A", "B", "C", "D"),
                                n_chains = 4, n_iter = 6000, burn_in = 2000,
                                seed = 42) {
  lambda_exp <- c(lambda1, lambda1 / 2)
  P_switch   <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  
  sim_exp <- simulate_vjm(v = v_true, lambda = lambda_exp, P = P_switch,
                          sigma = sigma_true, delta_t = delta_t,
                          T_total = T_total, seed = seed)
  
  prior_exp <- make_prior(n_states = 2, v_prior_sd = 5000,
                          log_lambda_mean = 0, log_lambda_sd = 3,
                          log_sigma_mean = log(sigma_true), log_sigma_sd = 1)
  
  cat(sprintf("\n========== lambda = (%.1f, %.1f) ==========\n",
              lambda_exp[1], lambda_exp[2]))
  # Switching summary for the regenerated dataset. switch_summary (Module 1)
  # returns both the expected rate under the stationary distribution and the
  # realised counts, including how many intervals exceed each truncation order.
  sw <- switch_summary(sim_exp)
  cat(sprintf("Switches per interval: %.2f expected, %.2f observed   (N = %d)\n",
              sw$expected, sw$per_interval, sw$n_intervals))
  cat(sprintf("Intervals with 0 / 1 / 2+ / 3+ switches: %d / %d / %d / %d\n",
              sw$n_0, sw$n_1, sw$n_2plus, sw$n_3plus))
  
  out <- list()
  for (L in likelihoods) {
    m <- run_mcmc(increments = sim_exp$increments, delta_t = delta_t,
                  prior = prior_exp, n_chains = n_chains, n_iter = n_iter,
                  burn_in = burn_in, seed = 123, likelihood = L,
                  verbose = FALSE)
    s <- posterior_summary(m)$natural_scale
    cat(sprintf("  %s: v1=%7.0f v2=%8.0f l1=%5.2f l2=%5.2f sigma=%5.1f | max R-hat=%.3f\n",
                L, s$mean[1], s$mean[2], s$mean[3], s$mean[4], s$mean[5],
                max(gelman_rubin(m))))
    out[[L]] <- m
  }
  cat(sprintf("  truth: v1=%7.0f v2=%8.0f l1=%5.2f l2=%5.2f sigma=%5.1f\n",
              v_true[1], v_true[2], lambda_exp[1], lambda_exp[2], sigma_true))
  invisible(out)
}


# SECTION 9: Experiment — sigma held at its true value

# All four likelihoods overestimate sigma (about 59 against a true 50). One
# hypothesis is that this inflated noise estimate makes the model insensitive
# to particle position, which in turn degrades the velocities and the rates.
# This experiment tests that directly by removing sigma from the parameter
# vector and holding it at the truth, then comparing against the same sampler
# with sigma free.

# Note that pinning sigma with a very tight prior is NOT a substitute: doing so
# destabilised the sampler in testing. `fix_sigma` removes the parameter.

run_fixed_sigma_experiment <- function(likelihood = "D", n_chains = 4,
                                       n_iter = 6000, burn_in = 2000) {
  cat(sprintf("\n========== sigma fixed vs free (Option %s) ==========\n",
              likelihood))
  for (fx in list(NULL, sigma_true)) {
    m <- run_mcmc(increments = sim$increments, delta_t = delta_t,
                  prior = prior, n_chains = n_chains, n_iter = n_iter,
                  burn_in = burn_in, seed = 123, likelihood = likelihood,
                  fix_sigma = fx, verbose = FALSE)
    s   <- posterior_summary(m)$natural_scale
    sg  <- if (is.null(fx)) sprintf("%5.1f", s$mean[5]) else sprintf("%5.1f*", fx)
    cat(sprintf("  sigma %-5s: v1=%7.0f v2=%8.0f l1=%5.2f l2=%5.2f sigma=%s | max R-hat=%.3f\n",
                if (is.null(fx)) "free" else "fixed",
                s$mean[1], s$mean[2], s$mean[3], s$mean[4], sg,
                max(gelman_rubin(m))))
  }
  cat(sprintf("  truth      : v1=%7.0f v2=%8.0f l1=%5.2f l2=%5.2f sigma=%5.1f\n",
              v_true[1], v_true[2], lambda_true[1], lambda_true[2], sigma_true))
  cat("  * held fixed, not estimated\n")
}


# Both experiments are defined above but not run automatically, since each
# takes several minutes. Run them explicitly when needed:
#run_rate_experiment(1)            # baseline, 0.3 switches/interval
#run_rate_experiment(4)          # faster switching (1.2 switches/interval)
#run_fixed_sigma_experiment("D") # sigma held at the truth

# SECTION 10: Experiment — thinning the observation series

# An alternative route to straining the truncated likelihoods, suggested as a
# counterpart to the rate experiment in Section 8: instead of switching faster,
# observe less often. Keeping every k-th observation multiplies the effective
# interval by k, so more switches fall between consecutive observations without
# changing the underlying process at all.

# Note the thinning is applied to the POSITIONS, which are then differenced.
# Thinning the increments themselves would be wrong: an increment is already a
# difference of two positions, so dropping every other increment discards
# displacement rather than lengthening the interval. Differencing the retained
# positions gives y_k - y_0, y_2k - y_k, ..., which is exactly what would have
# been recorded had the data been collected at delta_t * k.

# The noise structure is unchanged: each retained increment still differences
# two independently noisy positions, so its measurement variance remains
# 2 * sigma^2 and the emission models in Modules 2, 2B, 2C and 2D still apply.

# One caveat when interpreting the results: thinning changes two things at once.
# It lengthens the interval (the effect of interest) but also reduces the number
# of increments from N to roughly N / k, so some of any degradation is simply
# having less data. The rate experiment in Section 8 isolates the switching
# effect more cleanly and should be treated as the primary evidence, with this
# experiment as corroboration.

thin_observations <- function(sim, k) {
  idx <- seq(1, length(sim$positions), by = k)
  list(
    increments = diff(sim$positions[idx]),
    delta_t    = k * sim$params$delta_t,
    n          = length(idx) - 1
  )
}

run_thinning_experiment <- function(k = 2, likelihoods = c("A", "B", "C", "D"),
                                    n_chains = 4, n_iter = 6000,
                                    burn_in = 2000) {
  thinned <- thin_observations(sim, k)
  
  cat(sprintf("\n========== keeping every %d%s observation ==========\n",
              k, switch(as.character(k), "1" = "st", "2" = "nd", "3" = "rd", "th")))
  cat(sprintf("Effective delta_t: %.2f   (N = %d, was %d at delta_t = %.2f)\n",
              thinned$delta_t, thinned$n, length(sim$increments), delta_t))
  # Switching summary for the thinned series. The expected rate uses the
  # stationary distribution (proportional to the OPPOSITE rates) and the longer
  # effective interval; the observed counts are obtained by re-binning the exact
  # Gillespie jump times into the retained observation times.
  pi_stat  <- rev(lambda_true) / sum(lambda_true)
  exp_rate <- sum(pi_stat * lambda_true) * thinned$delta_t
  
  idx    <- seq(1, length(sim$positions), by = k)
  jt     <- sim$gillespie$jump_times
  jt     <- jt[jt > 0 & jt < sim$params$T_total]
  # [AI-GENERATED] Assignment of exact switch times to observation intervals
  # using findInterval and tabulate. See GenAI statement.
  counts <- tabulate(findInterval(jt, sim$times[idx], rightmost.closed = TRUE),
                     nbins = thinned$n)
  
  cat(sprintf("Switches per interval: %.2f expected, %.2f observed\n",
              exp_rate, length(jt) / thinned$n))
  cat(sprintf("Intervals with 0 / 1 / 2+ / 3+ switches: %d / %d / %d / %d\n",
              sum(counts == 0), sum(counts == 1),
              sum(counts >= 2), sum(counts >= 3)))
  
  for (L in likelihoods) {
    m <- run_mcmc(increments = thinned$increments, delta_t = thinned$delta_t,
                  prior = prior, n_chains = n_chains, n_iter = n_iter,
                  burn_in = burn_in, seed = 123, likelihood = L,
                  verbose = FALSE)
    s <- posterior_summary(m)$natural_scale
    cat(sprintf("  %s: v1=%7.0f v2=%8.0f l1=%5.2f l2=%5.2f sigma=%6.1f | max R-hat=%.3f\n",
                L, s$mean[1], s$mean[2], s$mean[3], s$mean[4], s$mean[5],
                max(gelman_rubin(m))))
  }
  cat(sprintf("  truth: v1=%7.0f v2=%8.0f l1=%5.2f l2=%5.2f sigma=%6.1f\n",
              v_true[1], v_true[2], lambda_true[1], lambda_true[2], sigma_true))
}

#run_thinning_experiment(2)

# SECTION 11: Experiment — larger sample size

# The rate bias at N = 200 survives every sampler, prior and likelihood
# diagnostic, which leaves sample size as the remaining explanation: the
# baseline dataset contains only about 50 switching events from which to
# estimate two rate parameters. This experiment lengthens the observation
# window, holding all model parameters and the sampling interval fixed, so the
# only thing that changes is how much information the data carry.

# Two tuning arguments matter at this sample size and are reported with the
# results, since they are sampler settings rather than model choices:
#   v_floor    the proposal SD floor for velocities. The default of 100 is set
#              for the N = 200 posterior; at N = 1600 the posterior SD of v1 is
#              roughly 15, so a floor of 100 prevents adaptation from shrinking
#              the proposal and collapses the acceptance rate.
#   Sigma_init the initial proposal covariance. The default diagonal is also
#              scaled for N = 200 (velocity SD 400), which is far wider than the
#              N = 1600 posterior, so burn-in is largely spent shrinking it.
#              Supplying a covariance matched to the expected posterior width
#              leaves more of the burn-in for adaptation to learn the shape.

run_sample_size_experiment <- function(T_exp = 480, v_floor = 5,
                                       Sigma_init = NULL,
                                       likelihoods = c("A", "B", "C", "D"),
                                       n_chains = 4, n_iter = 6000,
                                       burn_in = 2000, seed = 42) {
  P_switch <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  
  sim_exp <- simulate_vjm(v = v_true, lambda = lambda_true, P = P_switch,
                          sigma = sigma_true, delta_t = delta_t,
                          T_total = T_exp, seed = seed)
  
  sw <- switch_summary(sim_exp)
  cat(sprintf("\n========== T = %.0f  (N = %d) ==========\n",
              T_exp, sw$n_intervals))
  cat(sprintf("Switching events: %d   (%.2f per interval)\n",
              sw$n_switches, sw$per_interval))
  cat(sprintf("Velocity proposal floor: %.0f   Sigma_init: %s\n",
              v_floor, if (is.null(Sigma_init)) "default" else "supplied"))
  
  out <- list()
  for (L in likelihoods) {
    m <- run_mcmc(increments = sim_exp$increments, delta_t = delta_t,
                  prior = prior, n_chains = n_chains, n_iter = n_iter,
                  burn_in = burn_in, seed = 123, likelihood = L,
                  v_floor = v_floor, Sigma_init = Sigma_init,
                  verbose = FALSE)
    s <- posterior_summary(m)$natural_scale
    cat(sprintf("  %s: v1=%7.0f v2=%8.0f l1=%5.2f l2=%5.2f sigma=%5.1f | max R-hat=%.3f\n",
                L, s$mean[1], s$mean[2], s$mean[3], s$mean[4], s$mean[5],
                max(gelman_rubin(m))))
    out[[L]] <- m
  }
  cat(sprintf("  truth: v1=%7.0f v2=%8.0f l1=%5.2f l2=%5.2f sigma=%5.1f\n",
              v_true[1], v_true[2], lambda_true[1], lambda_true[2], sigma_true))
  invisible(out)
}

# The N = 1600 results reported in the dissertation were produced with:
#   run_sample_size_experiment(likelihoods = "A", n_chains = 6, n_iter = 15000, burn_in = 5000)
#   run_sample_size_experiment(likelihoods = "B", n_chains = 6, n_iter = 15000, burn_in = 5000)
#   run_sample_size_experiment(likelihoods = "C", n_chains = 6, n_iter = 15000, burn_in = 5000)
#   run_sample_size_experiment(likelihoods = "D", n_chains = 6, n_iter = 15000, burn_in = 5000)

# SECTION 12: Diagnostics reported in the results

# Three checks used to establish that the parameter bias at N = 200 is a
# property of the posterior rather than an artefact of the sampler, the
# starting point or the prior. Each is reported in the results; they are
# collected here so the reported numbers are reproducible.

# 12.1 Per-chain posterior means

# Where R-hat is inflated it is essential to distinguish a chain trapped in a
# secondary mode from chains that agree but have not fully equilibrated along a
# correlated direction. The maximum R-hat cannot make that distinction; the
# per-chain means can.

report_chains <- function(mcmc_out, label = "") {
  n    <- mcmc_out$prior$n_states
  post <- (mcmc_out$burn_in + 1):mcmc_out$n_iter
  
  cat(sprintf("\n--- Per-chain posterior means %s ---\n", label))
  cat(sprintf("%6s %8s %9s %7s %7s %8s\n",
              "chain", "v1", "v2", "lambda1", "lambda2", "sigma"))
  for (ch in seq_along(mcmc_out$chains)) {
    x <- mcmc_out$chains[[ch]]$chain[post, , drop = FALSE]
    cat(sprintf("%6d %8.0f %9.0f %7.2f %7.2f %8.1f\n", ch,
                mean(x[, 1]), mean(x[, 2]),
                mean(exp(x[, n + 1])), mean(exp(x[, n + 2])),
                if (is.null(mcmc_out$fix_sigma)) mean(exp(x[, 2 * n + 1]))
                else mcmc_out$fix_sigma))
  }
  cat(sprintf("%6s %8.0f %9.0f %7.2f %7.2f %8.1f\n", "truth",
              v_true[1], v_true[2], lambda_true[1], lambda_true[2], sigma_true))
  cat("\nPer-parameter R-hat:\n")
  print(round(gelman_rubin(mcmc_out), 3))
  invisible(NULL)
}

# 12.2 Start-at-truth

# All chains are initialised at the true parameters with negligible dispersion.
# If the chains remain there, an estimate away from the truth obtained from
# dispersed starts indicates incomplete convergence. If they drift away and
# settle elsewhere, that location is the posterior mode and the discrepancy is
# a property of the posterior rather than of the starting point.

run_start_at_truth <- function(likelihood = "D", n_chains = 4,
                               n_iter = 6000, burn_in = 2000) {
  theta_truth <- c(v_true, log(lambda_true), log(sigma_true))
  tiny_disp   <- c(5, 5, 0.02, 0.02, 0.02)
  
  m <- run_mcmc(increments = sim$increments, delta_t = delta_t, prior = prior,
                n_chains = n_chains, n_iter = n_iter, burn_in = burn_in,
                seed = 123, theta_centre = theta_truth, start_disp = tiny_disp,
                likelihood = likelihood, verbose = FALSE)
  
  s <- posterior_summary(m)$natural_scale
  cat(sprintf("\n--- Start-at-truth, Option %s ---\n", likelihood))
  cat(sprintf("  v1=%7.0f v2=%8.0f l1=%5.2f l2=%5.2f sigma=%5.1f | max R-hat=%.3f\n",
              s$mean[1], s$mean[2], s$mean[3], s$mean[4], s$mean[5],
              max(gelman_rubin(m))))
  cat(sprintf("  truth: v1=%7.0f v2=%8.0f l1=%5.2f l2=%5.2f sigma=%5.1f\n",
              v_true[1], v_true[2], lambda_true[1], lambda_true[2], sigma_true))
  invisible(m)
}

# 12.3 Prior sensitivity

# The rate prior is widened and then re-centred between the two true rates. If
# the posterior estimate is unchanged, the prior is not responsible for the
# bias and the likelihood is. Re-centring is the sharper test, since it moves
# the prior toward the truth and would pull the estimate with it if the prior
# were the cause.

run_prior_sensitivity <- function(likelihood = "D", n_chains = 4,
                                  n_iter = 6000, burn_in = 2000) {
  # The re-centred prior is placed at the geometric mean of the two true rates,
  # sqrt(lambda_1 * lambda_2) = 0.71, which is the value a prior centred between
  # them would take. It moves the prior toward the truth, so an estimate driven
  # by the prior would move with it.
  lambda_gm <- sqrt(prod(lambda_true))
  settings <- list(
    list(label = "default    (mean 0, sd 2)",  mean = 0,              sd = 2),
    list(label = "widened    (mean 0, sd 4)",  mean = 0,              sd = 4),
    list(label = "re-centred (mean log 0.71)", mean = log(lambda_gm), sd = 2)
  )
  
  cat(sprintf("\n--- Prior sensitivity on the rates, Option %s ---\n", likelihood))
  out <- list()
  for (cfg in settings) {
    pr <- make_prior(n_states = 2, log_lambda_mean = cfg$mean,
                     log_lambda_sd = cfg$sd,
                     log_sigma_mean = log(sigma_true), log_sigma_sd = 1)
    m  <- run_mcmc(increments = sim$increments, delta_t = delta_t, prior = pr,
                   n_chains = n_chains, n_iter = n_iter, burn_in = burn_in,
                   seed = 123, likelihood = likelihood, verbose = FALSE)
    s  <- posterior_summary(m)$natural_scale
    cat(sprintf("  %-27s lambda2 = %.2f | max R-hat = %.3f\n",
                cfg$label, s$mean[4], max(gelman_rubin(m))))
    out[[cfg$label]] <- m
  }
  cat(sprintf("  %-27s lambda2 = %.2f\n", "truth", lambda_true[2]))
  invisible(out)
}

# 12.4 Joint likelihood surface

# The joint surface over v1 and log sigma reported in Section 4.4 was computed
# for Option A alone. Because the Discussion appeals to its shape when
# explaining the behaviour of the other likelihoods, the surface is computed
# here for any requested option, so that the comparison can be made rather than
# assumed. `sim_obj` allows the surface to be evaluated on any of the simulated
# datasets, with the reference parameters taken from that dataset, so that the
# surface can be inspected at the switching rate where a given behaviour was
# observed rather than only at the baseline.

plot_joint_surface <- function(likelihood = "D", sim_obj = sim, n_grid = 60,
                               v1_range = c(1600, 2400),
                               logsig_range = c(3.5, 5.0),
                               file = NULL, drops = NULL) {
  ll_fun <- switch(likelihood,
                   "A" = log_likelihood_from_theta,
                   "B" = log_likelihood_from_theta_B,
                   "C" = log_likelihood_from_theta_C,
                   "D" = log_likelihood_from_theta_D,
                   stop("likelihood must be 'A', 'B', 'C' or 'D'"))
  
  pr_true    <- sim_obj$params
  theta_true <- c(pr_true$v, log(pr_true$lambda), log(pr_true$sigma))
  dt_obj     <- pr_true$delta_t
  
  v1_grid  <- seq(v1_range[1], v1_range[2], length.out = n_grid)
  sig_grid <- seq(logsig_range[1], logsig_range[2], length.out = n_grid)
  ll_2d    <- matrix(NA_real_, n_grid, n_grid)
  for (i in seq_along(v1_grid)) for (j in seq_along(sig_grid)) {
    th <- theta_true; th[1] <- v1_grid[i]; th[5] <- sig_grid[j]
    ll_2d[i, j] <- ll_fun(th, sim_obj$increments, dt_obj)
  }
  
  # Work in log-likelihood relative to the grid maximum. Absolute values differ
  # by hundreds of units between options, so relative values put every surface
  # on one scale and make the contour spacing interpretable as evidence.
  rel <- ll_2d - max(ll_2d, na.rm = TRUE)
  # Evenly spaced levels so the colour key is read linearly. Everything more
  # than `floor_drop` below the maximum is clipped to the lowest band; the
  # region of interest is the neighbourhood of the maximum, not the far tail.
  if (is.null(drops)) drops <- 40
  lev <- seq(-drops, 0, length.out = 11)
  rel_c <- pmax(rel, min(lev))
  
  idx  <- which(ll_2d == max(ll_2d, na.rm = TRUE), arr.ind = TRUE)[1, ]
  v1_h <- v1_grid[idx[1]]; sg_h <- sig_grid[idx[2]]
  
  if (!is.null(file)) png(file, width = 1150, height = 800, res = 130)
  op <- par(no.readonly = TRUE); on.exit(par(op))
  
  filled.contour(v1_grid, sig_grid, rel_c,
                 levels = lev,
                 col    = hcl.colors(length(lev) - 1, "Viridis", rev = TRUE),
                 xlab = expression(v[1]), ylab = expression(log(sigma)),
                 main = bquote("Joint log-likelihood surface: Option" ~ .(likelihood)),
                 key.title = title(main = "drop from\nmaximum", cex.main = 0.75, font.main = 1),
                 plot.axes = {
                   axis(1); axis(2)
                   contour(v1_grid, sig_grid, rel_c, levels = lev, add = TRUE,
                           col = "white", lwd = 0.7, drawlabels = FALSE)
                   # 2 log units marks the conventional support boundary; draw it heavier
                   contour(v1_grid, sig_grid, rel_c, levels = lev, add = TRUE,
                           col = "white", lwd = 0.7, drawlabels = FALSE)
                   abline(v = theta_true[1], col = "grey20", lty = 2, lwd = 1.2)
                   abline(h = theta_true[5], col = "grey20", lty = 2, lwd = 1.2)
                   points(theta_true[1], theta_true[5], pch = 4,  col = "grey10", cex = 1.6, lwd = 2.5)
                   points(v1_h, sg_h,                   pch = 21, col = "grey10", bg = "white",
                          cex = 1.5, lwd = 2)
                 })
  if (!is.null(file)) invisible(dev.off())
  
  cat(sprintf("Option %s: maximum at v1 = %.0f, sigma = %.1f; truth is %.2f log units below\n",
              likelihood, v1_h, exp(sg_h),
              max(ll_2d, na.rm = TRUE) - ll_fun(theta_true, sim_obj$increments, dt_obj)))
  invisible(list(v1 = v1_grid, log_sigma = sig_grid, ll = ll_2d, rel = rel))
}
#plot_joint_surface("A", file = "figures/fig8_joint_surface_A.png")
#plot_joint_surface("D", file = "figures/fig9_joint_surface_D.png")

#P <- matrix(c(0,1,1,0), 2, byrow = TRUE)
#sim8 <- simulate_vjm(v = v_true, lambda = c(8,4), P = P, sigma = sigma_true,
#                     delta_t = delta_t, T_total = T_total, seed = 42)
#plot_joint_surface("D", sim_obj = sim8, file = "figures/fig10_surface_D_lambda8.png")