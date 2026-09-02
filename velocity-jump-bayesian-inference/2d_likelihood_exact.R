# MODULE 2D: Exact-Switching Likelihood (Ceccarelli et al. P_infinity)

# GENERATIVE AI ATTRIBUTION
#   Passages developed with generative AI assistance are marked inline with
#   [AI-GENERATED] or [AI-ASSISTED] tags. [AI-GENERATED] marks code
#   substantially produced by AI; [AI-ASSISTED] marks code written by me and
#   subsequently debugged, corrected or optimised with AI assistance.
#   Unmarked code is my own. Full details are in the GenAI statement.

# PURPOSE:
#   Compute the log-likelihood of an observed increment track under an exact
#   treatment of within-interval switching, with no truncation of the switch
#   count. Each observation interval may contain any number of state switches,
#   and we integrate over the unknown switch times in every case.

#     Option A (Module 2)  = up-to-zero-switch   (m = 0)
#     Option B (Module 2b) = up-to-one-switch    (m = 1)
#     Option C (Module 2c) = up-to-two-switch    (m = 2)
#     Option D (this file) = exact               (m -> infinity)

#   Options A/B/C are successive truncations of Ceccarelli, Browning & Baker's
#   (2025, Bull. Math. Biol. 87:57) infinite sum over switch counts (their
#   Eq. 3). Option D is the limit of that sum, obtained in closed form via the
#   occupation-time distribution of the hidden state rather than by summing
#   switch-count terms. It combines the exact continuous-time dynamics of
#   Siekmann et al. (2011) with the velocity-jump emission of Ceccarelli et al.
#   (2025), and serves as the exact reference against which the truncated
#   Options A/B/C are judged.

# OCCUPATION-TIME GEOMETRY:
#   Over an interval of length delta_t the exact displacement depends on the
#   path only through T1, the total time spent in the start state s1:
#     delta_x = v_{s1} T1 + v_{s2} (delta_t - T1),
#   where s2 is the other state. Summing over all switch counts and switch
#   times is therefore equivalent to marginalising over the law of T1. For a
#   two-state CTMC started in s1, with l1 = lambda_{s1} (rate of leaving s1)
#   and l2 = lambda_{s2}, that law is
#     - a point mass at T1 = delta_t (no switch), weight exp(-l1 delta_t),
#       giving delta_x = v_{s1} delta_t and end state s1;
#     - a continuous density on u = T1 in (0, delta_t),
#         g(u) = exp(-l1 u - l2 (delta_t - u)) *
#                [ l1 I0(z) + sqrt(l1 l2 u / (delta_t - u)) I1(z) ],
#         z = 2 sqrt(l1 l2 u (delta_t - u)),
#       where I0, I1 are modified Bessel functions of the first kind.
#   (Verified by simulation, by exact normalisation, and by its one-switch
#   truncation reproducing the Module 2B density.)

# END-STATE ROUTING (important for the recursion):
#   The end state equals the start state s1 after an even number of switches
#   and the other state s2 after an odd number. The continuous density splits
#   accordingly:
#     end = s1 (even, >= 2 switches):  sqrt(l1 l2 u / w) I1(z)   (w = delta_t-u)
#     end = s2 (odd,  >= 1 switch)  :  l1 I0(z)
#   The even part carries little mass, since returning to s1 requires at least
#   two switches; the two parts sum to g(u) above. (Split verified against
#   end-state-resolved simulation.)

# EVALUATION:
#   The occupation-time density is closed-form, but the convolution with the
#   N(0, 2 sigma^2) measurement noise is evaluated by composite Simpson
#   quadrature over u in [0, delta_t] (mapped to delta_x). This is the same
#   scheme used in Module 2B and is accurate to the log-likelihood at the grid
#   sizes used here (513 points; convergence checked against grid refinement).

# INTERFACE:
#   log_likelihood_from_theta_D(theta, increments, delta_t, n_states = 2),
#   identical signature to Modules 2, 2B and 2C, so Module 3 can select it with
#   a single argument.


source("1_simulation.r")    # build_Q_matrix


# SECTION 1: Occupation-time emission densities, split by end state

# occ_even_odd returns the continuous occupation-time density on (0, delta_t),
# split by end state. l_start is the rate of leaving the interval's start
# state, l_other the rate of leaving the other state. "even" is the mass that
# ends back in the start state (>= 2 switches); "odd" is the mass that ends in
# the other state (>= 1 switch).
occ_even_odd <- function(u, l_start, l_other, delta_t) {
  ok   <- u > 0 & u < delta_t
  even <- numeric(length(u))
  odd  <- numeric(length(u))
  uu   <- u[ok]
  ww   <- delta_t - uu
  z    <- 2 * sqrt(l_start * l_other * uu * ww)
  base <- exp(-l_start * uu - l_other * ww)
  even[ok] <- base * sqrt(l_start * l_other * uu / ww) * besselI(z, 1)
  odd[ok]  <- base * l_start * besselI(z, 0)
  list(even = even, odd = odd)
}

# compute_emissions_D precomputes the emission tensor E[i, s1, s2], the density
# of increment i given start state s1 and end state s2, integrated exactly over
# occupation time. As in Module 2B the emissions depend only on the increment
# value and parameters, not on the forward variable, so they are computed once
# up front in vectorised form. The zero-switch atom contributes to E[, s1, s1];
# the continuous part is split between E[, s1, s1] (even) and E[, s1, s2] (odd).
compute_emissions_D <- function(increments, v, lambda, delta_t, sigma,
                                n_grid = 513) {
  n   <- length(v)
  N   <- length(increments)
  sdn <- sqrt(2) * sigma
  if (n_grid %% 2 == 0) n_grid <- n_grid + 1
  
  u  <- seq(0, delta_t, length.out = n_grid)
  # [AI-GENERATED] Composite Simpson weights by strided vector assignment.
  w  <- rep(2, n_grid); w[seq(2, n_grid - 1, by = 2)] <- 4; w[1] <- 1; w[n_grid] <- 1
  hu <- delta_t / (n_grid - 1)
  
  # [AI-GENERATED] Three-dimensional emission tensor and the partial array
  # indexing used to accumulate the no-switch atom and the two parity
  # components into the appropriate slices. See GenAI statement.
  E <- array(0, dim = c(N, n, n))
  for (s1 in 1:n) {
    s2 <- if (s1 == 1) 2 else 1
    l1 <- lambda[s1]; l2 <- lambda[s2]
    v1 <- v[s1];      v2 <- v[s2]
    
    # Zero-switch atom: no switch, end state s1, delta_x = v1 delta_t.
    atom_w <- exp(-l1 * delta_t)
    E[, s1, s1] <- E[, s1, s1] + atom_w * dnorm(increments, v1 * delta_t, sdn)
    
    # Continuous part, split by end state and convolved with Gaussian noise.
    eo <- occ_even_odd(u, l1, l2, delta_t)
    dx <- v1 * u + v2 * (delta_t - u)          # map occupation time to delta_x
    coef_even <- (hu / 3) * w * eo$even        # end = s1
    coef_odd  <- (hu / 3) * w * eo$odd         # end = s2
    # [AI-ASSISTED] begin: parity routing of the occupation-time density to end
    # states, with vectorised convolution. The first implementation transposed
    # the two Bessel terms, assigning the even and odd components to the wrong
    # end states. The error was identified by comparison against an
    # end-state-resolved Monte-Carlo simulation and corrected with AI assistance.
    # The form below is the corrected version. See GenAI statement.
    PHI <- dnorm(outer(increments, dx, FUN = "-"), mean = 0, sd = sdn)
    E[, s1, s1] <- E[, s1, s1] + as.vector(PHI %*% coef_even)
    E[, s1, s2] <- E[, s1, s2] + as.vector(PHI %*% coef_odd)
    # [AI-ASSISTED] end
  }
  E
}


# SECTION 2: Exact forward recursion

# The emission tensor E[i, s1, s2] already carries the exact joint law of
# (start state, end state, increment) over an interval, integrated over all
# switch counts and switch times. The recursion therefore routes mass from
# start state s1 to end state s2 through E directly, with no separate
# transition matrix and no switch-count truncation. This is the exact analogue
# of the exp(Q delta_t) routing in Module 2, generalised to the velocity-jump
# emission.
forward_algorithm_D <- function(increments, v, lambda, delta_t, sigma) {
  n <- length(v)
  N <- length(increments)
  
  # For a 2-state chain the only switch from each state is to the other, so
  # P = [[0,1],[1,0]]. For n > 2 we assume uniform switching to the other
  # states, matching the convention in Modules 2B and 2C.
  if (n == 2) {
    P_switch <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  } else {
    P_switch <- matrix(1, n, n); diag(P_switch) <- 0
    P_switch <- P_switch / rowSums(P_switch)
  }
  Q  <- build_Q_matrix(lambda, P_switch)   # from Module 1
  pi <- stationary_dist_D(Q)
  
  E <- compute_emissions_D(increments, v, lambda, delta_t, sigma)
  
  alpha   <- pi                            # interval 1 starts from stationary
  log_lik <- 0
  for (i in seq_len(N)) {
    new <- numeric(n)
    for (s1 in 1:n) {
      if (alpha[s1] <= 0) next
      for (s2 in 1:n) {
        new[s2] <- new[s2] + alpha[s1] * E[i, s1, s2]
      }
    }
    c_i <- sum(new)
    if (!is.finite(c_i) || c_i <= 0) return(-Inf)   # guard
    log_lik <- log_lik + log(c_i)
    alpha   <- new / c_i                            # renormalise (filtering)
  }
  log_lik
}

# Stationary distribution of Q. For two states the balance equations have the
# closed form pi = (lambda_2, lambda_1) / (lambda_1 + lambda_2), which is used
# directly: solving the system numerically instead fails when a proposed rate
# is large enough that the augmented matrix becomes numerically rank-deficient,
# which the sampler can reach under a wide rate prior. The general solve is
# retained for n > 2, guarded so that a failure returns the uniform
# distribution rather than aborting the run.
stationary_dist_D <- function(Q) {
  n <- nrow(Q)
  
  if (n == 2) {
    lambda <- -diag(Q)
    if (!all(is.finite(lambda)) || sum(lambda) <= 0) return(rep(1 / n, n))
    pi <- c(lambda[2], lambda[1])
    return(pi / sum(pi))
  }
  
  A  <- rbind(t(Q), rep(1, n))
  b  <- c(rep(0, n), 1)
  # [AI-GENERATED] Stationary distribution by qr.solve on the augmented balance
  # equations, with a fallback to a uniform distribution if the system is
  # numerically singular.
  pi <- tryCatch(as.vector(qr.solve(A, b)), error = function(e) rep(1 / n, n))
  pi <- pmax(pi, 0)
  pi / sum(pi)
}


# SECTION 3: theta wrapper

# theta = c(v_1, ..., v_n, log_lambda_1, ..., log_lambda_n, log_sigma),
# the same convention as Modules 2, 2B and 2C.
log_likelihood_from_theta_D <- function(theta, increments, delta_t,
                                        n_states = 2) {
  n <- n_states
  v      <- theta[1:n]
  lambda <- exp(theta[(n + 1):(2 * n)])
  sigma  <- exp(theta[2 * n + 1])
  
  if (any(!is.finite(c(v, lambda, sigma))) || any(lambda <= 0) || sigma <= 0) {
    return(-Inf)
  }
  forward_algorithm_D(increments, v, lambda, delta_t, sigma)
}


# SECTION 4: Validation (D against A/B/C at the truth)

# Confirms the exact likelihood is finite at the true parameters and that the
# emission integrates to one over the increment for each start state. Run only
# when this file is executed directly, not when sourced by Module 3.
# [AI-GENERATED] Guard controlling execution when the module is sourced.
if (sys.nframe() == 0) {
  source("2a_likelihood_zero_switch.r")   # Option A
  source("2b_likelihood_one_switch.r")    # Option B
  source("2c_likelihood_two_switch.r")    # Option C
  
  theta_true <- c(2000, -1500, log(1), log(0.5), log(50))
  
  ll_A <- log_likelihood_from_theta(theta_true,   sim$increments, delta_t)
  ll_B <- log_likelihood_from_theta_B(theta_true, sim$increments, delta_t)
  ll_C <- log_likelihood_from_theta_C(theta_true, sim$increments, delta_t)
  ll_D <- log_likelihood_from_theta_D(theta_true, sim$increments, delta_t)
  
  cat("=== MODULE 2D: LIKELIHOOD VALIDATION ===\n")
  
  # Log-likelihood and cost per evaluation at the true parameters. Option D is
  # not cheaper than the one-switch approximation, but costs about the same as
  # the two-switch approximation while being exact rather than truncated, so
  # the hierarchy crosses over at m = 2: a three-switch option would cost more
  # than D and still be approximate. D's cost is set by the emission quadrature
  # grid rather than by a switch count, so it does not grow with switching rate.
  lls <- c(A = ll_A, B = ll_B, C = ll_C, D = ll_D)
  fns <- list(A = log_likelihood_from_theta,   B = log_likelihood_from_theta_B,
              C = log_likelihood_from_theta_C, D = log_likelihood_from_theta_D)
  
  cat(sprintf("%-8s %14s %12s\n", "option", "log-lik", "ms/eval"))
  for (nm in names(fns)) {
    t0 <- Sys.time()
    for (i in 1:20) fns[[nm]](theta_true, sim$increments, delta_t)
    ms <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / 20 * 1000
    cat(sprintf("%-8s %14.2f %12.2f\n", nm, lls[[nm]], ms))
  }
  # Independent validation of the occupation-time emission.

  # The derivation is checked against three independent references: direct
  # Monte-Carlo simulation of CTMC paths, exact normalisation, and agreement
  # with the separately validated one-switch density of Module 2B. None of
  # these reuses the derivation being tested.
  
  l1 <- exp(theta_true[3]); l2 <- exp(theta_true[4])
  
  # Monte-Carlo reference: simulate paths directly and record the occupation
  # time in the start state together with the end state. This uses only the
  # exponential holding times, not the Bessel density.
  mc_occupation <- function(start_state, l1, l2, dt, n = 4e5, seed = 7) {
    set.seed(seed)
    lam <- c(l1, l2)
    T1 <- numeric(n); endst <- integer(n)
    for (i in seq_len(n)) {
      t <- 0; s <- start_state; time_in_start <- 0
      repeat {
        w <- rexp(1, lam[s])
        if (t + w >= dt) {
          if (s == start_state) time_in_start <- time_in_start + (dt - t)
          break
        }
        if (s == start_state) time_in_start <- time_in_start + w
        t <- t + w
        s <- if (s == 1) 2 else 1
      }
      T1[i] <- time_in_start; endst[i] <- s
    }
    list(T1 = T1, end = endst)
  }
  
  mc <- mc_occupation(1, l1, l2, delta_t)
  
  # Check 1: the no-switch atom. P(no switch) = exp(-l1 * delta_t).
  atom_analytic <- exp(-l1 * delta_t)
  atom_mc       <- mean(mc$T1 >= delta_t - 1e-9)
  
  # Checks 2 and 3: the continuous density, resolved by end state, against a
  # histogram of the simulated occupation times.
  brks <- seq(0, delta_t, length.out = 41)
  mids <- (brks[-1] + brks[-length(brks)]) / 2
  cont <- mc$T1 < delta_t - 1e-9
  eo   <- occ_even_odd(mids, l1, l2, delta_t)
  
  err <- numeric(2)
  for (es in 1:2) {
    sel <- cont & mc$end == es
    h   <- hist(mc$T1[sel], breaks = brks, plot = FALSE)
    emp <- h$density * mean(sel)
    ana <- if (es == 1) eo$even else eo$odd
    err[es] <- mean(abs(emp - ana))
  }
  
  # Check 4: the even and odd parts must sum to the total density g(u).
  uu  <- mids; ww <- delta_t - uu
  z   <- 2 * sqrt(l1 * l2 * uu * ww)
  tot <- exp(-l1 * uu - l2 * ww) *
    (l1 * besselI(z, 0) + sqrt(l1 * l2 * uu / ww) * besselI(z, 1))
  parity_err <- max(abs((eo$even + eo$odd) - tot))
  
  cat("\n--- Validation of the occupation-time emission ---\n")
  cat(sprintf("No-switch atom, analytic vs Monte Carlo : %.5f vs %.5f\n",
              atom_analytic, atom_mc))
  cat(sprintf("Occupation density vs MC, end = start   : mean abs error %.4f\n",
              err[1]))
  cat(sprintf("Occupation density vs MC, end = other   : mean abs error %.4f\n",
              err[2]))
  cat(sprintf("Even and odd parts sum to total density : max error %.3e\n",
              parity_err))
  
  # Check 5: emission should integrate to one over the increment for each start state.
  gy <- seq(-5000, 5000, length.out = 6001)
  E  <- compute_emissions_D(gy, theta_true[1:2],
                            exp(theta_true[3:4]), delta_t, exp(theta_true[5]))
  for (s1 in 1:2) {
    total <- sum(E[, s1, 1] + E[, s1, 2]) * (gy[2] - gy[1])
    cat(sprintf("Emission integral, start s%d = %.4f  (should be ~1)\n",
                s1, total))
  }
  
  # Check 6: truncating the switch count at one must reproduce the Module 2B
  # one-switch density. The exactly-one-switch term of the occupation density
  # is l1 * exp(-l1 u) * exp(-l2 (delta_t - u)), so the k <= 1 emission is the
  # no-switch atom plus that term convolved with the noise. Agreement is
  # checked as the quadrature grid is refined; a monotone decrease confirms
  # the residual is discretisation error and not an error in the derivation.
  # [AI-ASSISTED] Degenerate case guard for equal rates, using a tolerance
  # rather than exact equality. See GenAI statement.
  Z_one <- function(l1, l2, dt) {
    if (abs(l1 - l2) < 1e-12) l1 * exp(-l1 * dt) * dt
    else l1 * exp(-l2 * dt) * (1 - exp(-(l1 - l2) * dt)) / (l1 - l2)
  }
  emission_k1 <- function(delta_y, v, lambda, dt, sigma, n_grid) {
    v1 <- v[1]; v2 <- v[2]; ls1 <- lambda[1]; ls2 <- lambda[2]
    sdn <- sqrt(2) * sigma
    if (n_grid %% 2 == 0) n_grid <- n_grid + 1
    u  <- seq(0, dt, length.out = n_grid)
    gu <- ifelse(u > 0 & u < dt, ls1 * exp(-ls1 * u) * exp(-ls2 * (dt - u)), 0)
    dx <- v1 * u + v2 * (dt - u)
    # [AI-GENERATED] Composite Simpson weights by strided vector assignment.
    w  <- rep(2, n_grid); w[seq(2, n_grid - 1, 2)] <- 4
    w[1] <- 1; w[n_grid] <- 1
    coef <- (dt / (n_grid - 1) / 3) * w * gu
    exp(-ls1 * dt) * dnorm(delta_y, v1 * dt, sdn) +
      # [AI-GENERATED] Vectorised convolution. See GenAI statement.
      as.vector(dnorm(outer(delta_y, dx, "-"), 0, sdn) %*% coef)
  }
  
  dy_test <- seq(-3000, 3500, length.out = 14)
  v_t     <- theta_true[1:2]; lam_t <- c(l1, l2); sg_t <- exp(theta_true[5])
  ref_2B  <- exp(-l1 * delta_t) * dnorm(dy_test, v_t[1] * delta_t, sqrt(2) * sg_t) +
    Z_one(l1, l2, delta_t) *
    one_switch_density(dy_test, 1, 2, v_t, lam_t, delta_t, sg_t,
                       n_grid = 1025)
  
  cat("One-switch truncation vs Module 2B:\n")
  for (ng in c(201, 513, 1025, 4097)) {
    k1 <- emission_k1(dy_test, v_t, lam_t, delta_t, sg_t, ng)
    cat(sprintf("   n_grid = %5d : max relative difference %.4f\n",
                ng, max(abs(k1 - ref_2B) / pmax(ref_2B, 1e-12))))
  }
  
  # Quadrature-grid sensitivity.

  # Option D's cost is set by the emission quadrature grid rather than by a
  # switch count. Refining the grid must converge the log-likelihood; the rate
  # of convergence shows how coarse a grid is tolerable, and therefore how far
  # the cost can be reduced without material loss of accuracy.
  
  ll_D_grid <- function(theta, n_grid) {
    v      <- theta[1:2]
    lambda <- exp(theta[3:4])
    sigma  <- exp(theta[5])
    P      <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
    Q      <- build_Q_matrix(lambda, P)
    alpha  <- stationary_dist_D(Q)
    E      <- compute_emissions_D(sim$increments, v, lambda, delta_t, sigma,
                                  n_grid = n_grid)
    log_lik <- 0
    for (i in seq_along(sim$increments)) {
      new <- numeric(2)
      for (s1 in 1:2) for (s2 in 1:2) new[s2] <- new[s2] + alpha[s1] * E[i, s1, s2]
      c_i <- sum(new)
      if (!is.finite(c_i) || c_i <= 0) return(-Inf)
      log_lik <- log_lik + log(c_i)
      alpha   <- new / c_i
    }
    log_lik
  }
  
  cat("\n--- Option D: quadrature-grid sensitivity ---\n")
  ref <- ll_D_grid(theta_true, 2049)
  cat(sprintf("%8s %14s %14s %10s\n",
              "n_grid", "log-lik", "diff vs ref", "ms/eval"))
  for (ng in c(65, 129, 257, 513, 1025)) {
    ll <- ll_D_grid(theta_true, ng)
    t0 <- Sys.time()
    for (i in 1:5) ll_D_grid(theta_true, ng)
    ms <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / 5 * 1000
    cat(sprintf("%8d %14.4f %14.2e %10.2f\n", ng, ll, abs(ll - ref), ms))
  }
}