# MODULE 2B: Up-to-One-Switch Likelihood (Ceccarelli et al. P_1)

# GENERATIVE AI ATTRIBUTION
#   Passages developed with generative AI assistance are marked inline with
#   [AI-GENERATED] or [AI-ASSISTED] tags. [AI-GENERATED] marks code
#   substantially produced by AI; [AI-ASSISTED] marks code written by me and
#   subsequently debugged, corrected or optimised with AI assistance.
#   Unmarked code is my own. Full details are in the GenAI statement.

# PURPOSE:
#   Compute the log-likelihood of an observed increment track under the
#   "up-to-one-switch" approximation of Ceccarelli et al. (2025),
#   "Approximate solutions of a general stochastic velocity-jump model
#   subject to discrete-time noisy observations" (Bull. Math. Biol. 87:57, 2025).

#   Option B is the up-to-one-switch (m = 1) member of the up-to-m-switch
#   hierarchy. The three options compared in this project are:
#     Option A (Module 2)   up-to-zero-switch: each interval is treated as pure
#                           single-state motion; biased when switching is
#                           frequent relative to delta_t.
#     Option B (Module 2b)  up-to-one-switch: each interval may contain zero or
#                           one switch, integrating over the unknown switch time.
#     Option C (Module 2c)  up-to-two-switch.

#   Both A and B are members of Ceccarelli's up-to-m-switch hierarchy
#   (their Equation 5): A is m = 0, B is m = 1.

# WHY OPTION A IS BIASED:
#   When a switch occurs mid-interval, the true displacement is
#     delta_x = v_{s1} * tau1 + v_{s2} * (delta_t - tau1),
#   i.e. somewhere between v_{s1}*delta_t and v_{s2}*delta_t, depending on the
#   random switch time tau1. Option A has no way to represent such an
#   increment except as noise, so it inflates sigma and shrinks the velocity
#   gap. Option B represents these intervals correctly by integrating over
#   tau1.

# DATA MODEL (shared with Modules 1, 2):
#   delta_y_j = delta_x_j + delta_eps_j,  delta_eps_j ~ N(0, 2*sigma^2).
#   The variance is 2*sigma^2 because the increment differences two
#   independent N(0, sigma^2) measurement errors.

# JOINT vs MARGINAL:
#   Ceccarelli Section 3 gives the marginal single-increment PDF P_1(delta_y).
#   Ceccarelli Section 4 promotes it to the JOINT PDF over a whole track via a
#   forward recursion that carries the filtered start-of-interval state
#   distribution P(S_1^j | delta_y_{1:j-1}) from one interval to the next.
#   We implement the JOINT version here, because (a) it is the correct
#   replacement for the Option-A forward algorithm, and (b) Module 3 
#   already expects a single track log-likelihood.

#   The recursion exploits the fact that the state at the END of interval j is
#   the state at the START of interval j+1 (Ceccarelli Figure 8). Under the
#   up-to-one-switch approximation:
#     - if zero switches occur in interval j, the end state equals the start
#       state s1;
#     - if one switch occurs, the end state is the second state s2.
#   We assume, as Ceccarelli do, that the measurement-noise increments of
#   successive intervals are independent (their Equation S4 assumption).

# INTERFACE:
#   log_likelihood_from_theta_B(theta, increments, delta_t, n_states = 2)
#   with theta = c(v_1..v_n, log_lambda_1..log_lambda_n, log_sigma),
#   identical to Module 2's log_likelihood_from_theta(). Module 3 can switch
#   between A and B by changing only the function name it calls.


source("1_simulation.r")   # simulation + build_Q_matrix


# SECTION 1: CTMC building blocks

# Stationary distribution pi of the CTMC, solving pi Q = 0.
# For a 2-state chain with rates lambda1, lambda2 this is
#   pi = (lambda2, lambda1) / (lambda1 + lambda2),
# but we compute it generally from Q so the code extends to n states.
stationary_dist <- function(Q) {
  n <- nrow(Q)
  
  # Closed form for the common 2-state case: pi proportional to (lambda2,
  # lambda1) where lambda_s = -Q[s,s]. Robust and exact.
  if (n == 2) {
    lambda <- -diag(Q)
    pi <- c(lambda[2], lambda[1])
    return(pi / sum(pi))
  }
  
  # General case: pi is the left null vector of Q (solves pi Q = 0), i.e. the
  # eigenvector of Q^T for eigenvalue 0.
  # [AI-GENERATED] Stationary distribution by eigendecomposition of t(Q).
  e   <- eigen(t(Q))
  idx <- which.min(abs(e$values))
  pi  <- Re(e$vectors[, idx])
  pi  <- pi / sum(pi)
  pi[pi < 0] <- 0
  pi <- pi / sum(pi)
  return(pi)
}

# Embedded DTMC switching probabilities p_{su} = q_{su} / lambda_s for u != s,
# p_{ss} = 0. For a 2-state chain every off-diagonal p is 1.
embedded_P <- function(Q) {
  n <- nrow(Q)
  lambda <- -diag(Q)
  P <- matrix(0, n, n)
  for (s in 1:n) {
    if (lambda[s] > 0) {
      for (u in 1:n) {
        if (u != s) P[s, u] <- Q[s, u] / lambda[s]
      }
    }
  }
  return(P)
}


# SECTION 2: One-switch increment density f_tilde_{s1,s2}(delta_y)

# Density of a noisy increment given exactly one switch from state s1 to s2
# within the interval. Derived (and numerically validated against simulation)
# as follows:

#   Exact increment given switch time tau1:
#       delta_x = v_{s2}*delta_t + (v_{s1} - v_{s2}) * tau1,
#   which is linear and monotonic in tau1, so delta_x ranges over
#   [min(v1 dt, v2 dt), max(v1 dt, v2 dt)].

#   Density of tau1 given one switch (s1 then s2), from the joint exponential
#   holding times with tau2 > delta_t - tau1, normalised over (0, delta_t]:
#       p(tau1) propto lambda_{s1} exp(-lambda_{s1} tau1)
#                                   * exp(-lambda_{s2} (delta_t - tau1)).

#   Change of variables tau1 -> delta_x (Jacobian 1/|v1 - v2|) gives the exact
#   increment density g_tilde(delta_x). Convolving with N(0, 2 sigma^2) gives
#   the noisy density f_tilde(delta_y).

# We evaluate f_tilde by 1-D Simpson quadrature of the convolution over the
# compact support [min(v1,v2)*dt, max(v1,v2)*dt]. The tau1 normaliser Z is
# computed in closed form. This is accurate (log-likelihood converged at 129
# grid points, verified against adaptive quadrature and direct simulation) and
# fast enough for MCMC. Ceccarelli's SI S4 gives a fully closed-form f_tilde
# that could replace the convolution quadrature for additional speed.

one_switch_density <- function(delta_y, s1, s2, v, lambda, delta_t, sigma,
                               n_grid = 129) {
  v1 <- v[s1]; v2 <- v[s2]
  l1 <- lambda[s1]; l2 <- lambda[s2]
  sd_noise <- sqrt(2) * sigma
  
  # If the two velocities are equal, a switch produces no kink: the increment
  # is just N(v1 delta_t, 2 sigma^2) (Ceccarelli special case v_s1 = v_s2).
  # [AI-ASSISTED] Degenerate case guard. Written by me and corrected with AI
  # assistance to use all.equal rather than exact equality, since floating-point
  # velocities proposed by the sampler are never exactly equal. See GenAI statement.
  if (isTRUE(all.equal(v1, v2))) {
    return(dnorm(delta_y, mean = v1 * delta_t, sd = sd_noise))
  }
  
  # Normalising constant Z for the tau1 density over (0, delta_t].
  # Z = INT_0^dt l1 exp(-l1 t1) exp(-l2 (dt - t1)) dt1, which has the closed
  # form below (avoids an inner adaptive integration on every call, which is
  # important because this function runs millions of times inside MCMC):
  #   l1 != l2 :  l1 exp(-l2 dt) (1 - exp(-(l1-l2) dt)) / (l1 - l2)
  #   l1 == l2 :  l1 exp(-l1 dt) dt
  # [AI-GENERATED] begin: closed-form normalising constant for the switch-time
  # density, including the separate branch required when the two rates are
  # equal. See GenAI statement.
  if (abs(l1 - l2) < 1e-12) {
    Z <- l1 * exp(-l1 * delta_t) * delta_t
  } else {
    Z <- l1 * exp(-l2 * delta_t) * (1 - exp(-(l1 - l2) * delta_t)) / (l1 - l2)
  }
  # [AI-GENERATED] end
  tau1_kernel <- function(t1) l1 * exp(-l1 * t1) * exp(-l2 * (delta_t - t1))
  
  # [AI-GENERATED] begin: choice of integration domain and quadrature grid.
  # The exact-increment density is compactly supported, so the convolution is
  # integrated only over that support; the Gaussian tail outside it is captured
  # by evaluating the noise density at grid points inside it. Integrating beyond
  # the support would straddle the discontinuity at the edges and lose accuracy.
  # The composite Simpson weights are built by strided vector assignment.
  # See GenAI statement.
  # Support of the exact-increment density g_tilde: [lo, hi].
  # g_tilde is compactly supported on [lo, hi] and zero outside, so we
  # integrate the convolution only over [lo, hi]; the Gaussian tail for
  # delta_y outside the support is captured automatically by evaluating
  # phi(delta_y - dx) at dx in [lo, hi]. Integrating beyond [lo, hi] is both
  # unnecessary and harmful, because g_tilde is discontinuous at the edges and
  # a grid straddling the jump loses accuracy.
  lo <- min(v1 * delta_t, v2 * delta_t)
  hi <- max(v1 * delta_t, v2 * delta_t)
  
  # Composite Simpson's rule over [lo, hi] (n_grid must be odd; 129 points
  # gives a log-likelihood converged to >4 d.p., verified against both
  # adaptive quadrature and direct simulation).
  if (n_grid %% 2 == 0) n_grid <- n_grid + 1
  dx_grid <- seq(lo, hi, length.out = n_grid)
  h       <- (hi - lo) / (n_grid - 1)
  
  # g_tilde on the grid (interior is smooth; endpoints are the finite jump)
  t1 <- (dx_grid - v2 * delta_t) / (v1 - v2)
  gd <- tau1_kernel(t1) / Z / abs(v1 - v2)
  
  # Simpson weights: 1, 4, 2, 4, ..., 4, 1
  w        <- rep(2, n_grid)
  w[seq(2, n_grid - 1, by = 2)] <- 4
  w[1] <- 1; w[n_grid] <- 1
  # [AI-GENERATED] end
  
  # Convolve with Gaussian noise for every delta_y at once.
  # f_tilde(dy_i) = sum_k (h/3) w_k gd_k phi(dy_i - dx_k).
  # Writing coef_k = (h/3) w_k gd_k, this is a matrix-vector product
  #   f = PHI %*% coef,  where PHI[i,k] = phi(dy_i - dx_k).
  # Vectorising over delta_y this way is ~15-20x faster than sapply-looping,
  # with an identical result.
  # [AI-GENERATED] begin: vectorised convolution. The integral of the emission
  # density against the measurement-noise density is expressed as an
  # outer-difference matrix multiplied by a coefficient vector, rather than
  # looping over increments. See GenAI statement.
  coef <- (h / 3) * w * gd
  # outer difference dy_i - dx_k, then Gaussian density, then matrix multiply
  PHI  <- dnorm(outer(delta_y, dx_grid, FUN = "-"), mean = 0, sd = sd_noise)
  as.vector(PHI %*% coef)
  # [AI-GENERATED] end
}

# SECTION 3: Joint forward recursion (Ceccarelli Section 4)

# Carries alpha[s] = P(S_1^j = s, delta_y_{1:j-1}) — the unnormalised
# probability that interval j starts in state s, given all increments so far.
# At each interval we split the mass by start state, apply the up-to-one-
# switch densities, route the mass to the correct end state (= next start
# state), accumulate the normaliser into the log-likelihood, and renormalise.

forward_algorithm_B <- function(increments, v, lambda, delta_t, sigma) {
  n  <- length(v)
  N  <- length(increments)
  
  # build_Q_matrix (Module 1) requires the embedded switching matrix P as well
  # as the rates. For a 2-state chain the only possible switch from each state
  # is to the other, so P = [[0,1],[1,0]]. For n > 2 with no further
  # information we assume uniform switching to the other states; this can be
  # generalised later when the inference includes the p_{su}.
  if (n == 2) {
    P_switch <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  } else {
    P_switch <- matrix(1, n, n)
    diag(P_switch) <- 0
    P_switch <- P_switch / rowSums(P_switch)
  }
  
  Q  <- build_Q_matrix(lambda, P_switch)   # from Module 1
  pi <- stationary_dist(Q)
  P  <- embedded_P(Q)
  
  # Per-start-state switch probabilities for one interval
  pW0   <- exp(-lambda * delta_t)     # P(no switch | start state s)
  pWge1 <- 1 - pW0                    # P(at least one switch | start state s)
  
  # ---- PRECOMPUTE EMISSION TENSORS (vectorised over ALL increments) ----
  # The emission densities depend only on the increment value and parameters,
  # not on the forward variable alpha, so they can be computed once up front in
  # vectorised form rather than re-evaluated inside the recursion. This is the
  # key optimisation: it turns ~N*n*n separate density calls (each rebuilding a
  # quadrature grid) into n*n vectorised calls, giving a large speed-up with an
  # identical log-likelihood (verified to machine precision).
  #
  #   E0[i, s]      = zero-switch density of increment i given start state s
  #   E1[i, s1, s2] = one-switch density of increment i given s1 -> s2
  sd_noise <- sqrt(2) * sigma
  
  E0 <- matrix(0, nrow = N, ncol = n)
  for (s in 1:n) {
    E0[, s] <- dnorm(increments, mean = v[s] * delta_t, sd = sd_noise)
  }
  
  # [AI-GENERATED] Three-dimensional emission array indexed by increment,
  # start state and end state.
  E1 <- array(0, dim = c(N, n, n))
  for (s1 in 1:n) {
    for (s2 in 1:n) {
      if (s1 == s2) next
      # one_switch_density is vectorised: pass the whole increment vector
      E1[, s1, s2] <- one_switch_density(increments, s1, s2,
                                         v, lambda, delta_t, sigma)
    }
  }
  
  # ---- CHEAP RECURSION over precomputed emissions ----
  alpha   <- pi                       # interval 1 starts from stationary dist
  log_lik <- 0
  
  for (i in seq_len(N)) {
    new <- numeric(n)
    
    for (s1 in 1:n) {
      if (alpha[s1] <= 0) next
      
      # Zero-switch: stays in s1; end state = s1
      new[s1] <- new[s1] + alpha[s1] * pW0[s1] * E0[i, s1]
      
      # One-switch: s1 -> s2; end state = s2
      for (s2 in 1:n) {
        if (s2 == s1) next
        new[s2] <- new[s2] + alpha[s1] * pWge1[s1] * P[s1, s2] * E1[i, s1, s2]
      }
    }
    
    c_i <- sum(new)
    if (!is.finite(c_i) || c_i <= 0) return(-Inf)   # guard
    
    log_lik <- log_lik + log(c_i)
    alpha   <- new / c_i                            # renormalise (filtering)
  }
  
  return(log_lik)
}


# SECTION 4: theta wrapper (same convention as Module 2)

# theta = c(v_1, ..., v_n, log_lambda_1, ..., log_lambda_n, log_sigma)
log_likelihood_from_theta_B <- function(theta, increments, delta_t,
                                        n_states = 2) {
  n <- n_states
  v      <- theta[1:n]
  lambda <- exp(theta[(n + 1):(2 * n)])
  sigma  <- exp(theta[2 * n + 1])
  
  if (any(!is.finite(c(v, lambda, sigma))) || any(lambda <= 0) || sigma <= 0) {
    return(-Inf)
  }
  
  forward_algorithm_B(increments, v, lambda, delta_t, sigma)
}


# SECTION 5: Validation against the Option-A likelihood

# Compares Option B with Option A at the true parameters and at perturbations.
# The key result we expect: under Option B the true parameters should sit much
# closer to the peak of the likelihood, so a far larger fraction of random
# perturbations should decrease the log-likelihood (fraction-below-true should
# rise from ~0.4 toward ~1).

# [AI-GENERATED] Guard controlling execution when the module is sourced.
if (sys.nframe() == 0) {     # only run when sourced as a script
  
  cat("\n=== MODULE 2B: UP-TO-ONE-SWITCH LIKELIHOOD VALIDATION ===\n")
  
  source("2a_likelihood_zero_switch.r")   # brings in Option-A log_likelihood_from_theta
  
  theta_true <- c(v_true[1], v_true[2],
                  log(lambda_true[1]), log(lambda_true[2]),
                  log(sigma_true))
  
  ll_A <- log_likelihood_from_theta(theta_true, sim$increments, delta_t)
  ll_B <- log_likelihood_from_theta_B(theta_true, sim$increments, delta_t)
  
  cat(sprintf("\nLog-likelihood at TRUE theta:\n"))
  cat(sprintf("  Option A (zero-switch): %.2f\n", ll_A))
  cat(sprintf("  Option B (one-switch) : %.2f\n", ll_B))
  cat("  (B should be HIGHER: it explains switch intervals it used to mis-model)\n")
  
  # Perturbation test for Option B
  set.seed(1)
  perturb_sd <- c(200, 200, 0.3, 0.3, 0.3)
  n_perturb  <- 20
  ll_pB <- numeric(n_perturb)
  for (k in 1:n_perturb) {
    tp <- theta_true + rnorm(length(theta_true), sd = perturb_sd)
    ll_pB[k] <- log_likelihood_from_theta_B(tp, sim$increments, delta_t)
  }
  cat(sprintf("\nOption B perturbation test (%d perturbations):\n", n_perturb))
  cat(sprintf("  True LL  : %.2f\n", ll_B))
  cat(sprintf("  Mean pert: %.2f\n", mean(ll_pB)))
  cat(sprintf("  Max pert : %.2f\n", max(ll_pB)))
  cat(sprintf("  Fraction below true: %.2f  (want close to 1)\n",
              mean(ll_pB < ll_B)))
}