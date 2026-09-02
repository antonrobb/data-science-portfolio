# MODULE 2C: Up-to-Two-Switch Likelihood (Ceccarelli et al. P_2)

# GENERATIVE AI ATTRIBUTION
#   Passages developed with generative AI assistance are marked inline with
#   [AI-GENERATED] or [AI-ASSISTED] tags. [AI-GENERATED] marks code
#   substantially produced by AI; [AI-ASSISTED] marks code written by me and
#   subsequently debugged, corrected or optimised with AI assistance.
#   Unmarked code is my own. Full details are in the GenAI statement.

# PURPOSE:
#   Compute the log-likelihood of an observed increment track under the
#   "up-to-two-switch" approximation. It extends Option B by allowing
#   each observation interval to contain zero, one, or two state switches,
#   integrating over the unknown switch times in each case.

#     Option A (Module 2)  = up-to-zero-switch  (m = 0)
#     Option B (Module 2b) = up-to-one-switch   (m = 1)
#     Option C (this file) = up-to-two-switch   (m = 2)

#   These are successive truncations of Ceccarelli, Browning & Baker's
#   (2025, Bull. Math. Biol. 87:57) exact infinite sum over switch counts
#   (their Eq. 3). As m increases the
#   approximation improves; the gain from each extra order shrinks because
#   higher switch counts are exponentially less probable. With the present
#   parameters two-switch intervals occur only ~1-2% of the time, so Option C
#   contributes only a small correction relative to Option B.

# TWO-SWITCH GEOMETRY (the new piece):
#   In a 2-state model a two-switch interval goes s1 -> s2 -> s3 with each
#   consecutive pair distinct, which forces s3 = s1: the particle returns to
#   its starting state. With switch times 0 < tau1 < tau2 < delta_t the exact
#   displacement is
#     delta_x = v_{s1} tau1 + v_{s2} (tau2 - tau1) + v_{s1} (delta_t - tau2)
#             = v_{s1} delta_t + (v_{s2} - v_{s1}) * w,   w = tau2 - tau1.
#   So delta_x depends only on w, the total time spent in the other state s2.
#   Marginalising the joint exponential switch-time density over the start time
#   tau1 (with l3 = l1 because s3 = s1) gives the closed-form w-density
#     p(w) = l2 exp(-l2 w) exp(-l1 delta_t) l1 (delta_t - w),   w in (0, delta_t),
#   where l1 = lambda_{s1}, l2 = lambda_{s2}. (Verified against simulation.)
#   Changing variables w -> delta_x (Jacobian 1/|v2 - v1|) and convolving with
#   the N(0, 2 sigma^2) measurement noise yields the noisy two-switch density.

# END-STATE ROUTING (important for the recursion):
#   - 0 switches: end state = s1            (mass stays in s1)
#   - 1 switch  : end state = s2            (mass moves to s2)
#   - 2 switches: end state = s1            (mass returns to s1)
#   The forward recursion must route the two-switch contribution back to the
#   start state, unlike the one-switch contribution.

# SWITCH-COUNT PROBABILITIES (per start state s1, other state s2):
#   P0 = exp(-l1 delta_t)
#   P1 = INT_0^delta_t l1 e^{-l1 t1} e^{-l2 (delta_t - t1)} dt1            (= Z1)
#   P2 = INT_0^delta_t l2 e^{-l2 w} e^{-l1 delta_t} l1 (delta_t - w) dw
#   These do not sum to 1 (they omit >= 3 switches), exactly as in Ceccarelli's
#   truncated hierarchy; the missing mass here is < 0.4% and is intentionally
#   not renormalised, so Options A/B/C are directly comparable as successive
#   truncations.

# INTERFACE:
#   log_likelihood_from_theta_C(theta, increments, delta_t, n_states = 2),
#   identical signature to Modules 2 and 2B, so Module 3 can select it with a
#   single argument.


source("1_simulation.r")    # build_Q_matrix
source("2b_likelihood_one_switch.r")   # reuse stationary_dist, embedded_P,
# one_switch_density (and Option B for compare)


# SECTION 1: Closed-form switch-count probabilities

# Z1: the one-switch normaliser / probability (same closed form as Module 2B).
# [AI-GENERATED] begin: closed-form switch-count normalising constants,
# including the separate branches required when the two rates are equal.
# See GenAI statement.
Z1_prob <- function(l1, l2, delta_t) {
  if (abs(l1 - l2) < 1e-12) {
    l1 * exp(-l1 * delta_t) * delta_t
  } else {
    l1 * exp(-l2 * delta_t) * (1 - exp(-(l1 - l2) * delta_t)) / (l1 - l2)
  }
}

# P2: probability of exactly two switches s1 -> s2 -> s1 within the interval.
#   INT_0^dt l2 e^{-l2 w} e^{-l1 dt} l1 (dt - w) dw
# Closed form (verified by symbolic integration and against simulation):
#   P2 = l1 (dt l2 e^{dt l2} - e^{dt l2} + 1) e^{-dt (l1 + l2)} / l2
P2_prob <- function(l1, l2, delta_t) {
  e2 <- exp(delta_t * l2)
  l1 * (delta_t * l2 * e2 - e2 + 1) * exp(-delta_t * (l1 + l2)) / l2
}
# [AI-GENERATED] end

# SECTION 2: Two-switch increment density f2_{s1}(delta_y)

# Noisy density of an increment given exactly two switches starting (and
# ending) in state s1. delta_x ranges over [min(v1,v2)*dt, max(v1,v2)*dt] just
# like the one-switch case, but with a different (triangular-weighted) shape.
# Vectorised over delta_y by the same Simpson-convolution scheme as Module 2B.

two_switch_density <- function(delta_y, s1, v, lambda, delta_t, sigma,
                               n_grid = 129) {
  s2 <- if (s1 == 1) 2 else 1          # the other state (2-state model)
  v1 <- v[s1]; v2 <- v[s2]
  l1 <- lambda[s1]; l2 <- lambda[s2]
  sd_noise <- sqrt(2) * sigma
  
  # Equal-velocity degenerate case: a switch produces no kink.
  # [AI-ASSISTED] Degenerate case guard. Written by me and corrected with AI
  # assistance to use all.equal rather than exact equality, since floating-point
  # velocities proposed by the sampler are never exactly equal. See GenAI statement.
  if (isTRUE(all.equal(v1, v2))) {
    return(dnorm(delta_y, mean = v1 * delta_t, sd = sd_noise))
  }
  
  if (n_grid %% 2 == 0) n_grid <- n_grid + 1
  lo <- min(v1 * delta_t, v2 * delta_t)
  hi <- max(v1 * delta_t, v2 * delta_t)
  dx_grid <- seq(lo, hi, length.out = n_grid)
  h       <- (hi - lo) / (n_grid - 1)
  
  # w = time spent in state s2 = (dx - v1 dt) / (v2 - v1), must lie in (0, dt).
  # The unnormalised w-density is l2 e^{-l2 w} e^{-l1 dt} l1 (dt - w); it
  # integrates to P2 (the probability of exactly two switches), not to 1. To
  # use it as a conditional density f2(dy | W=2) -- matching how the one-switch
  # density is conditional on W=1 -- we divide by P2 so it integrates to 1. The
  # switch-count weight P(W>=2) is then applied separately in the recursion
  # (Ceccarelli Eq. 5).
  P2 <- P2_prob(l1, l2, delta_t)
  w  <- (dx_grid - v1 * delta_t) / (v2 - v1)
  gd <- ifelse(w > 0 & w < delta_t,
               l2 * exp(-l2 * w) * exp(-l1 * delta_t) * l1 * (delta_t - w) /
                 abs(v2 - v1) / P2,
               0)
  
  # Simpson weights
  # [AI-GENERATED] Vectorised convolution. The integral of the emission
  # density against the measurement-noise density is expressed as an
  # outer-difference matrix multiplied by a coefficient vector, rather than
  # looping over increments. See GenAI statement.
  wgt <- rep(2, n_grid); wgt[seq(2, n_grid - 1, by = 2)] <- 4
  wgt[1] <- 1; wgt[n_grid] <- 1
  coef <- (h / 3) * wgt * gd
  
  PHI <- dnorm(outer(delta_y, dx_grid, FUN = "-"), mean = 0, sd = sd_noise)
  as.vector(PHI %*% coef)
}


# SECTION 3: Joint forward recursion (up to two switches)

forward_algorithm_C <- function(increments, v, lambda, delta_t, sigma) {
  n <- length(v)
  N <- length(increments)
  
  if (n != 2) {
    stop("Module 2C (two-switch) is implemented for the 2-state model only. ",
         "For n > 2 the two-switch end-state routing has more cases.")
  }
  
  P_switch <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  Q  <- build_Q_matrix(lambda, P_switch)
  pi <- stationary_dist(Q)
  
  sd_noise <- sqrt(2) * sigma
  
  # Precompute switch-count probabilities per start state
  # Ceccarelli Eq. 5 with m = 2: the up-to-two-switch approximation is
  #   P0 * f0  +  P1 * f1  +  P(W>=2) * f2,
  # i.e. the lower orders (w = 0, 1) carry their exact probabilities and the
  # top order (w = 2) carries the tail probability P(W >= 2). This mirrors how
  # Option B (m = 1) weights its one-switch term by P(W >= 1), and is what
  # makes the hierarchy a consistent family of truncations.
  P0   <- exp(-lambda * delta_t)                   # exactly zero switches
  P1   <- numeric(n); Pge2 <- numeric(n)
  for (s1 in 1:n) {
    s2 <- if (s1 == 1) 2 else 1
    P1[s1]   <- Z1_prob(lambda[s1], lambda[s2], delta_t)   # exactly one
    Pge2[s1] <- 1 - P0[s1] - P1[s1]                        # two or more
  }
  
  # Precompute emission tensors (vectorised over all increments)
  # E0[i, s]      zero-switch density (end state s)
  # E1[i, s1, s2] one-switch density  (s1 -> s2, end state s2)
  # E2[i, s1]     two-switch density  (s1 -> s2 -> s1, end state s1)
  E0 <- matrix(0, N, n)
  for (s in 1:n) E0[, s] <- dnorm(increments, mean = v[s] * delta_t, sd = sd_noise)
  
  # [AI-GENERATED] Three-dimensional emission array indexed by increment,
  # start state and end state.
  E1 <- array(0, dim = c(N, n, n))
  for (s1 in 1:n) for (s2 in 1:n) {
    if (s1 == s2) next
    E1[, s1, s2] <- one_switch_density(increments, s1, s2,
                                       v, lambda, delta_t, sigma)
  }
  
  E2 <- matrix(0, N, n)
  for (s1 in 1:n) {
    E2[, s1] <- two_switch_density(increments, s1, v, lambda, delta_t, sigma)
  }
  
  # Cheap recursion over precomputed emissions
  alpha   <- pi
  log_lik <- 0
  
  for (i in seq_len(N)) {
    new <- numeric(n)
    for (s1 in 1:n) {
      if (alpha[s1] <= 0) next
      s2 <- if (s1 == 1) 2 else 1
      
      # 0 switches: end state s1
      new[s1] <- new[s1] + alpha[s1] * P0[s1] * E0[i, s1]
      # exactly 1 switch: end state s2
      new[s2] <- new[s2] + alpha[s1] * P1[s1] * E1[i, s1, s2]
      # 2 or more switches (approximated by the two-switch density): returns to s1
      new[s1] <- new[s1] + alpha[s1] * Pge2[s1] * E2[i, s1]
    }
    
    c_i <- sum(new)
    if (!is.finite(c_i) || c_i <= 0) return(-Inf)
    
    log_lik <- log_lik + log(c_i)
    alpha   <- new / c_i
  }
  
  log_lik
}


# SECTION 4: theta wrapper

log_likelihood_from_theta_C <- function(theta, increments, delta_t,
                                        n_states = 2) {
  n <- n_states
  v      <- theta[1:n]
  lambda <- exp(theta[(n + 1):(2 * n)])
  sigma  <- exp(theta[2 * n + 1])
  
  if (any(!is.finite(c(v, lambda, sigma))) || any(lambda <= 0) || sigma <= 0) {
    return(-Inf)
  }
  forward_algorithm_C(increments, v, lambda, delta_t, sigma)
}


# SECTION 5: Validation (A vs B vs C at the truth)

# [AI-GENERATED] Guard controlling execution when the module is sourced.
if (sys.nframe() == 0) {
  
  cat("\n=== MODULE 2C: UP-TO-TWO-SWITCH LIKELIHOOD VALIDATION ===\n")
  
  source("2a_likelihood_zero_switch.r")   # Option A
  
  theta_true <- c(v_true[1], v_true[2],
                  log(lambda_true[1]), log(lambda_true[2]),
                  log(sigma_true))
  
  ll_A <- log_likelihood_from_theta(theta_true,   sim$increments, delta_t)
  ll_B <- log_likelihood_from_theta_B(theta_true, sim$increments, delta_t)
  ll_C <- log_likelihood_from_theta_C(theta_true, sim$increments, delta_t)
  
  cat(sprintf("\nLog-likelihood at TRUE theta:\n"))
  cat(sprintf("  Option A (zero-switch): %.2f\n", ll_A))
  cat(sprintf("  Option B (one-switch) : %.2f\n", ll_B))
  cat(sprintf("  Option C (two-switch) : %.2f\n", ll_C))
  cat("  (Note: P_m is an APPROXIMATION, so the likelihood need not increase\n")
  cat("   monotonically in m at a fixed point. The relevant test is PDF\n")
  cat("   accuracy vs the true density and posterior recovery of sigma, not\n")
  cat("   the raw likelihood ordering. C should recover sigma slightly closer\n")
  cat("   to the truth than B.)\n")
  
  set.seed(1)
  perturb_sd <- c(200, 200, 0.3, 0.3, 0.3)
  n_perturb  <- 20
  ll_pC <- numeric(n_perturb)
  for (k in 1:n_perturb) {
    tp <- theta_true + rnorm(length(theta_true), sd = perturb_sd)
    ll_pC[k] <- log_likelihood_from_theta_C(tp, sim$increments, delta_t)
  }
  cat(sprintf("\nOption C perturbation test (%d perturbations):\n", n_perturb))
  cat(sprintf("  True LL  : %.2f\n", ll_C))
  cat(sprintf("  Mean pert: %.2f\n", mean(ll_pC)))
  cat(sprintf("  Max pert : %.2f\n", max(ll_pC)))
  cat(sprintf("  Fraction below true: %.2f\n", mean(ll_pC < ll_C)))
}