# MODULE 1: Simulation of a 2-State Velocity-Jump Model

# GENERATIVE AI ATTRIBUTION
#   Passages developed with generative AI assistance are marked inline with
#   [AI-GENERATED] or [AI-ASSISTED] tags. [AI-GENERATED] marks code
#   substantially produced by AI; [AI-ASSISTED] marks code written by me and
#   subsequently debugged, corrected or optimised with AI assistance.
#   Unmarked code is my own. Full details are in the GenAI statement.

# PURPOSE:
#   Simulate particle trajectories from a continuous-time Markov chain (CTMC)
#   with 2 hidden velocity states (Forward and Backward). The particle moves
#   at a constant velocity determined by its hidden state, switching between
#   states stochastically. Noisy position observations are recorded at fixed
#   time intervals delta_t.

# BACKGROUND:
#   This follows the velocity-jump model framework of Ceccarelli et al.
#   (2025, Bull. Math. Biol. 87:57) and the CTMC machinery of Siekmann et al.
#   (2011). The hidden state evolves as a CTMC with generator matrix Q.
#   The Gillespie algorithm is used to simulate the exact continuous-time
#   state sequence, which is then sampled at discrete intervals with
#   Gaussian measurement noise.

# STATE DEFINITIONS (2-state model):
#   State 1: Forward  (velocity v1 > 0)
#   State 2: Backward (velocity v2 < 0)

# PARAMETERS:
#   v      : numeric vector of length n, velocities for each state [units/time]
#   lambda : numeric vector of length n, total switching rates [1/time]
#   P      : n x n matrix, switching probabilities (rows sum to 1, diag = 0)
#   sigma  : scalar, standard deviation of Gaussian measurement noise
#   delta_t: scalar, fixed time between observations
#   T_total: scalar, total duration of the simulation

# OUTPUTS:
#   A list containing:
#     $times       : observation times t_0, t_1, ..., t_N
#     $positions   : noisy observed positions y_0, y_1, ..., y_N
#     $true_pos    : true (noiseless) positions at observation times
#     $true_states : hidden state at each observation time
#     $increments  : observed position increments delta_y_1, ..., delta_y_N
#     $gillespie   : full Gillespie output (state switches, exact times)


# SECTION 1: Build the Generator Matrix Q

# The generator matrix Q is the core mathematical object of a CTMC.
# For state i:
#   - Off-diagonal entry q_ij = lambda_i * P_ij  (rate of jumping to state j)
#   - Diagonal entry q_ii = -lambda_i             (total rate of leaving state i)

# This satisfies the conservative property: rows of Q sum to zero.
# Reference: Siekmann et al. (2011), Equation 2.

build_Q_matrix <- function(lambda, P) {
  # INPUT VALIDATION
  n <- length(lambda)
  stopifnot(
    "lambda must be a positive numeric vector" = is.numeric(lambda) && all(lambda > 0),
    "P must be a square matrix"               = is.matrix(P) && nrow(P) == ncol(P),
    "P must have n rows matching lambda"      = nrow(P) == n,
    "P diagonal must be zero"                 = all(diag(P) == 0),
    "P rows must sum to 1"                    = all(abs(rowSums(P) - 1) < 1e-10)
  )
  
  # Build off-diagonal entries: q_ij = lambda_i * p_ij
  Q <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        Q[i, j] <- lambda[i] * P[i, j]
      }
    }
  }
  
  # Set diagonal so rows sum to zero (conservative property)
  diag(Q) <- -rowSums(Q)
  
  return(Q)
}


# SECTION 2: Gillespie Algorithm for CTMC Simulation

# Simulates the exact continuous-time trajectory of the hidden state.

# The Gillespie algorithm works as follows:
#   1. In state i, the time until the next switch is Exponential(lambda_i)
#   2. When a switch occurs, jump to state j with probability P_ij
#   3. Repeat until total time exceeds T_total

# This gives the exact stochastic path of the hidden state — no approximation.
# Reference: Gillespie (1976), as cited in Siekmann et al. (2011).

simulate_gillespie <- function(Q, T_total, initial_state = 1) {
  n <- nrow(Q)
  
  # Switching rates: lambda_i = -Q_ii (total rate of leaving state i)
  lambda <- -diag(Q)
  
  # Storage for state sequence and jump times
  states     <- initial_state   # vector of states visited
  jump_times <- 0               # times at which state changes occur
  
  current_state <- initial_state
  current_time  <- 0
  
  while (current_time < T_total) {
    # Time until next switch: Exponential with rate lambda[current_state]
    wait_time <- rexp(1, rate = lambda[current_state])
    next_time <- current_time + wait_time
    
    # If next switch would exceed T_total, stop here
    if (next_time >= T_total) {
      jump_times <- c(jump_times, T_total)
      states     <- c(states, current_state)  # stays in current state until end
      break
    }
    
    # Choose which state to jump to using transition probabilities
    # Conditional jump probabilities: p_ij = Q_ij / lambda_i  (for j != i)
    jump_probs <- Q[current_state, ] / lambda[current_state]
    jump_probs[current_state] <- 0  # cannot stay in same state
    jump_probs <- pmax(jump_probs, 0)  # numerical safety: avoid tiny negatives
    
    next_state <- sample(1:n, size = 1, prob = jump_probs)
    
    # Record the jump
    current_time  <- next_time
    current_state <- next_state
    
    jump_times <- c(jump_times, current_time)
    states     <- c(states, current_state)
  }
  
  return(list(
    states     = states,      # state at each interval [initial, after jump 1, ...]
    jump_times = jump_times   # times of state changes [0, t1, t2, ..., T_total]
  ))
}


# SECTION 3: Compute True Position from Gillespie Output

# Given the exact Gillespie path, compute the particle's true position at any
# time t by integrating velocity over each state interval.

# Position at time t:
#   x(t) = x(0) + integral_0^t v(s(u)) du

# Since velocity is piecewise constant, this is a sum of (velocity * duration)
# over each Gillespie interval up to time t.

get_true_position <- function(gillespie, v, query_times, x0 = 0) {
  states     <- gillespie$states
  jump_times <- gillespie$jump_times
  
  n_intervals <- length(states) - 1  # number of complete intervals
  # Note: states[i] is active during [jump_times[i], jump_times[i+1])
  
  n_query     <- length(query_times)
  positions   <- numeric(n_query)
  true_states <- integer(n_query)
  
  for (q in 1:n_query) {
    t <- query_times[q]
    pos <- x0
    
    if (t == 0) {
      true_states[q] <- states[1]
      positions[q]   <- x0
      next
    }
    
    for (i in 1:n_intervals) {
      t_start <- jump_times[i]
      t_end   <- jump_times[i + 1]
      state_i <- states[i]
      
      if (t <= t_start) break  # query time is before this interval
      
      # Duration of this interval that overlaps with [0, t]
      effective_end <- min(t_end, t)
      duration      <- effective_end - t_start
      
      pos <- pos + v[state_i] * duration
      
      if (t <= t_end) {
        true_states[q] <- state_i
        break
      }
      
      # If we reach the last interval
      if (i == n_intervals) {
        true_states[q] <- states[i]
      }
    }
    
    positions[q] <- pos
  }
  
  return(list(positions = positions, states = true_states))
}


# SECTION 4: Main Simulation Function

# Ties everything together:
#   1. Build Q from parameters
#   2. Run Gillespie to get exact hidden path
#   3. Sample true positions at observation times
#   4. Add Gaussian noise to get observed positions

simulate_vjm <- function(v,
                         lambda,
                         P,
                         sigma,
                         delta_t,
                         T_total,
                         initial_state = 1,
                         x0            = 0,
                         seed          = NULL) {
  # Set random seed for reproducibility
  if (!is.null(seed)) set.seed(seed)
  
  # Step 1: Validate inputs
  n <- length(v)
  stopifnot(
    "v must be numeric"                  = is.numeric(v),
    "lengths of v and lambda must match" = length(lambda) == n,
    "P must have n rows and columns"     = all(dim(P) == c(n, n)),
    "sigma must be positive"             = sigma > 0,
    "delta_t must be positive"           = delta_t > 0,
    "T_total must be > delta_t"          = T_total > delta_t
  )
  
  # Step 2: Build Q matrix
  Q <- build_Q_matrix(lambda, P)
  
  # Step 3: Simulate hidden state path (Gillespie)
  gillespie <- simulate_gillespie(Q, T_total, initial_state)
  
  # Step 4: Define observation times
  # t_0 = 0, t_1 = delta_t, ..., t_N = N * delta_t
  obs_times <- seq(0, T_total, by = delta_t)
  N         <- length(obs_times) - 1  # number of increments
  
  # Step 5: Get true positions at observation times
  true_result  <- get_true_position(gillespie, v, obs_times, x0)
  true_pos     <- true_result$positions
  true_states  <- true_result$states
  
  # Step 6: Add Gaussian measurement noise
  noise        <- rnorm(length(obs_times), mean = 0, sd = sigma)
  observed_pos <- true_pos + noise
  
  # Step 7: Compute observed position increments
  # delta_y_i = y_i - y_{i-1}  for i = 1, ..., N
  increments <- diff(observed_pos)
  
  # Step 8: Return everything
  return(list(
    # Simulation parameters (store for reference)
    params = list(
      v       = v,
      lambda  = lambda,
      P       = P,
      Q       = Q,
      sigma   = sigma,
      delta_t = delta_t,
      T_total = T_total,
      n_states = n
    ),
    # Observation grid
    times        = obs_times,        # length N+1
    # Observed (noisy) data
    positions    = observed_pos,     # length N+1
    increments   = increments,       # length N  — THIS IS THE KEY DATA FOR INFERENCE
    # Ground truth (for validation)
    true_pos     = true_pos,         # length N+1
    true_states  = true_states,      # length N+1
    # Full Gillespie output
    gillespie    = gillespie
  ))
}

# Switching-event summary

# The number of switching events in a dataset is the quantity the rate
# parameters are estimated from, so it is worth reporting alongside N. It also
# determines how much the truncated likelihoods (Options A, B, C) are strained:
# the mean number of switches per observation interval indicates how often an
# interval contains more switches than the truncation allows.

# `n_switches` counts state changes in the exact Gillespie path. Because the
# jump times are known exactly, each switch can be assigned to the observation
# interval it falls in, giving the distribution of switches per interval. The
# intervals holding two or more switches are the ones Option B cannot represent,
# and those holding three or more are beyond Option C, so these counts say
# directly how much work the truncations are being asked to do.
switch_summary <- function(sim) {
  jt <- sim$gillespie$jump_times
  jt <- jt[jt > 0 & jt < sim$params$T_total]   # interior switches only
  N  <- length(sim$increments)
  
  # Assign each switch to its observation interval, then tabulate.
  # [AI-GENERATED] Assignment of exact switch times to observation intervals
  # using findInterval and tabulate. See GenAI statement.
  bin    <- findInterval(jt, sim$times, rightmost.closed = TRUE)
  counts <- tabulate(bin, nbins = N)
  
  list(
    n_switches   = length(jt),
    n_intervals  = N,
    per_interval = length(jt) / N,
    n_0          = sum(counts == 0),
    n_1          = sum(counts == 1),
    n_2plus      = sum(counts >= 2),
    n_3plus      = sum(counts >= 3)
  )
}

# SECTION 5: Visualisation Functions

# Plot 1: Full trajectory showing true path, noisy observations, hidden state

# The true path is piecewise linear, with a constant velocity between switches
# and a kink at each switch, so it is drawn as line segments taken directly from
# the exact Gillespie jump times rather than as points at the observation times.
# The noisy observations are drawn first and the truth superimposed, so that the
# true path is seen to lie in the middle of the observation scatter.

# At the parameters used for inference the measurement noise (sigma = 50) is
# invisible against a position range of tens of thousands, so `sigma_plot`
# allows the figure to be redrawn with an inflated noise level purely for
# illustration. The underlying true path is unchanged; only the observation
# noise is redrawn. Any figure produced this way is labelled as illustrative.
plot_trajectory <- function(sim, max_time = NULL, sigma_plot = NULL,
                            seed_plot = 1) {
  if (is.null(max_time)) max_time <- sim$params$T_total
  
  # Observation-level quantities
  idx         <- sim$times <= max_time
  times_plot  <- sim$times[idx]
  true_plot   <- sim$true_pos[idx]
  states_plot <- sim$true_states[idx]
  
  # Observations: either as simulated, or redrawn with an inflated sigma.
  if (is.null(sigma_plot)) {
    obs_plot   <- sim$positions[idx]
    sigma_used <- sim$params$sigma
    illustrative <- FALSE
  } else {
    set.seed(seed_plot)
    obs_plot   <- true_plot + rnorm(length(true_plot), 0, sigma_plot)
    sigma_used <- sigma_plot
    illustrative <- TRUE
  }
  
  state_cols  <- c("1" = "orange2", "2" = "dodgerblue3")
  state_names <- c("1" = "Forward", "2" = "Backward")
  n_states    <- sim$params$n_states
  
  # Exact true path from the Gillespie output: positions at the jump times,
  # joined by straight segments coloured by the state holding over each segment.
  gp <- sim$gillespie
  jt <- gp$jump_times
  js <- gp$states
  x0 <- sim$true_pos[1]
  tp <- get_true_position(gp, sim$params$v, jt, x0 = x0)$positions
  
  par(mfrow = c(2, 1), mar = c(5, 7, 3, 1))
  
  # Top panel: observations first, true path superimposed.
  yr     <- range(c(obs_plot, true_plot))
  yr_top <- yr[2] + 0.18 * diff(yr)
  plot(times_plot, obs_plot,
       type = "p", pch = 16, cex = 0.5, col = "grey55",
       xlab = "Time", ylab = "Position",
       main = "Simulated particle trajectory",
       sub  = if (illustrative)
         sprintf("Observation noise inflated to sigma = %.0f for visibility (inference uses sigma = %.0f)",
                 sigma_used, sim$params$sigma)
       else sprintf("sigma = %.0f", sigma_used),
       cex.sub = 0.75,
       ylim = c(yr[1], yr_top))
  
  for (i in seq_len(length(jt) - 1)) {
    if (jt[i] >= max_time) break
    t_end <- min(jt[i + 1], max_time)
    # linear interpolation if the segment is clipped by max_time
    x_end <- tp[i] + (tp[i + 1] - tp[i]) *
      (t_end - jt[i]) / (jt[i + 1] - jt[i])
    segments(jt[i], tp[i], t_end, x_end,
             col = state_cols[as.character(js[i])], lwd = 2)
  }
  
  legend("top", horiz = TRUE,
         legend = c("Observed (noisy)", paste("True:", state_names)),
         col    = c("grey55", state_cols),
         lty    = c(NA, 1, 1), pch = c(16, NA, NA),
         lwd    = c(NA, 2, 2), cex = 0.8, bty = "n")
  
  # Bottom panel: hidden state as an exact step function.
  plot(times_plot, states_plot,
       type = "s", col = "seagreen4", lwd = 2,
       xlab = "Time", ylab = "", main = "Hidden state sequence",
       yaxt = "n", ylim = c(0.5, n_states + 0.5))
  axis(2, at = 1:n_states,
       labels = paste0("State ", 1:n_states, " (", state_names, ")"),
       las = 1, cex.axis = 0.8)
  
  par(mfrow = c(1, 1))
}


# Plot 2: Distribution of observed increments
# Key for validating that the likelihood approximation is sensible
plot_increment_distribution <- function(sim) {
  inc <- sim$increments
  n_states <- sim$params$n_states
  v       <- sim$params$v
  sigma   <- sim$params$sigma
  delta_t <- sim$params$delta_t
  
  sd_increment <- sqrt(2) * sigma
  comp_means   <- v[1:n_states] * delta_t
  x_lo <- min(comp_means, inc) - 4 * sd_increment
  x_hi <- max(comp_means, inc) + 4 * sd_increment
  
  hist(inc,
       breaks = 50,
       freq   = FALSE,
       col    = "grey80",
       border = "white",
       xlab   = expression(Delta * y),
       main   = "Distribution of Observed Position Increments",
       sub    = paste0("N = ", length(inc), " increments"),
       xlim   = c(x_lo, x_hi),
       xaxt   = "n",
       ylim   = c(0, max(dnorm(0, mean = 0, sd = sqrt(2) * sigma)) * 1.2))
  
  axis(1, at = pretty(c(x_lo, x_hi), n = 10))
  
  # Overlay expected Gaussian components for each state
  # If no switching occurs in [t_{i-1}, t_i], then:
  #   delta_y ~ N(v_s * delta_t, 2 * sigma^2)
  x_seq <- seq(x_lo, x_hi, length.out = 500)
  cols   <- c("orange2", "dodgerblue3", "seagreen4", "orchid3")
  ltys   <- c(2, 4, 3, 6)
  
  # The OBSERVED increment differences two independent noisy positions, so
  # its noise has variance 2*sigma^2 (sd = sqrt(2)*sigma), matching the
  # emission model used in the likelihood (Modules 2, 2B, 2C). Using sd =
  # sigma here would draw the components too narrow to sit on the data.
  for (s in 1:n_states) {
    mean_s <- v[s] * delta_t
    lines(x_seq, dnorm(x_seq, mean = mean_s, sd = sd_increment),
          col = cols[s], lwd = 2, lty = ltys[s])
  }
  
  legend("topright",
         legend = paste0("State ", 1:n_states, ": N(", round(v * delta_t, 1), ", ", 2 * sigma^2, ")"),
         col    = cols[1:n_states],
         lwd    = 2, lty = ltys[1:n_states], bty = "n")
}

# Plot 3: Repeated realisations of the same model

# The velocity-jump model is stochastic, so a single parameter set does not
# determine a single trajectory. Simulating the same parameters several times
# shows how much the observed path can vary purely through the randomness of
# the switching times, which is the variability the inference has to work
# against. This is a methods figure: it makes clear why the parameters cannot
# be read directly off one track, and why a likelihood over the whole track is
# needed rather than a fit to any individual path.

# `n_reps` realisations are drawn using consecutive seeds, all with identical
# parameters. True paths are drawn as lines; if `show_obs` is TRUE the noisy
# observations of the first realisation are added as points, to show the
# measurement noise on top of the process noise.
plot_repeat_trajectories <- function(v, lambda, P, sigma, delta_t, T_total,
                                     n_reps = 5, seed0 = 100, show_obs = FALSE) {
  reps <- lapply(seq_len(n_reps), function(i) {
    simulate_vjm(v = v, lambda = lambda, P = P, sigma = sigma,
                 delta_t = delta_t, T_total = T_total,
                 initial_state = 1, x0 = 0, seed = seed0 + i)
  })
  
  cols <- c("orange2", "dodgerblue3", "seagreen4", "orchid3", "firebrick3")
  ltys <- c(1, 2, 3, 4, 5)
  
  yr <- range(unlist(lapply(reps, function(r) r$true_pos)))
  if (show_obs) yr <- range(c(yr, reps[[1]]$positions))
  yr[2] <- yr[2] + 0.20 * diff(yr)   # headroom so the legend clears the paths
  
  par(mfrow = c(1, 1), mar = c(4, 5, 3, 1))
  plot(NA,
       xlim = c(0, T_total),
       ylim = yr,
       xlab = "Time",
       ylab = "Position",
       main = sprintf("%d realisations of the same model", n_reps))
  abline(h = 0, col = "grey85")
  
  if (show_obs) {
    points(reps[[1]]$times, reps[[1]]$positions,
           pch = 16, cex = 0.4, col = "grey60")
  }
  
  for (i in seq_len(n_reps)) {
    lines(reps[[i]]$times, reps[[i]]$true_pos,
          col = cols[(i - 1) %% length(cols) + 1],
          lty = ltys[(i - 1) %% length(ltys) + 1], lwd = 2)
  }
  
  legend("top", horiz = TRUE,
         legend = seq_len(n_reps),
         title  = "Realisation",
         col    = cols[1:n_reps],
         lty    = ltys[1:n_reps],
         lwd    = 2, bty = "n", cex = 0.85)
  
  # Report how much the endpoint varies across realisations, as a compact
  # numerical companion to the visual spread.
  ends <- sapply(reps, function(r) tail(r$true_pos, 1))
  cat(sprintf("Final position across %d realisations: min %.0f, max %.0f, sd %.0f\n",
              n_reps, min(ends), max(ends), sd(ends)))
  
  invisible(reps)
}

# SECTION 6: Run the Simulation

# These parameters match the 2-state example from Ceccarelli et al. (2025)
# Figure 2, making it easy to compare results.

# True parameter values (theta_true)
# Following Ceccarelli et al. (2025), Figure 2B
v_true      <- c(2000, -1500)   # State 1: Forward, State 2: Backward [units/time]
lambda_true <- c(1.0,   0.5)    # Switching rates [1/time]
sigma_true  <- 50               # Measurement noise standard deviation

# Switching probability matrix P
# For 2 states: the only option is to switch to the other state, so P_12 = P_21 = 1
P_true <- matrix(c(0, 1,
                   1, 0),
                 nrow = 2, byrow = TRUE)

# Observation settings
delta_t <- 0.3    # time between observations (same as Ceccarelli et al.)
T_total <- 60     # total simulation time (gives N = 200 increments)

# Run simulation. Guarded by exists("sim") so that the simulation is generated
# and reported once per session: Modules 2, 2B, 2C, 2D and 3 each source this
# file, and without the guard the same dataset would be regenerated and printed
# on every source. The guard preserves the objects the downstream modules need,
# which a sys.nframe() guard would not.
if (!exists("sim")) {
  cat("Running simulation...\n")
  sim <- simulate_vjm(
    v             = v_true,
    lambda        = lambda_true,
    P             = P_true,
    sigma         = sigma_true,
    delta_t       = delta_t,
    T_total       = T_total,
    initial_state = 1,
    x0            = 0,
    seed          = 42          # set seed for reproducibility
  )
  
  # Print summary
  cat("\n=== SIMULATION SUMMARY ===\n")
  cat("Number of observations (N+1):    ", length(sim$times), "\n")
  cat("Number of increments (N):        ", length(sim$increments), "\n")
  sw <- switch_summary(sim)
  cat("Switching events:                ", sw$n_switches, "\n")
  cat("Mean switches per interval:      ", round(sw$per_interval, 3), "\n")
  cat("Intervals with 0 / 1 / 2+ / 3+:  ",
      sprintf("%d / %d / %d / %d", sw$n_0, sw$n_1, sw$n_2plus, sw$n_3plus), "\n")
  cat("True Q matrix:\n")
  print(round(sim$params$Q, 4))
  cat("\nState frequencies (hidden):\n")
  print(table(sim$true_states))
  cat("\nObserved increment summary:\n")
  print(summary(sim$increments))
}

# Plot results. Only when this file is run directly, so that sourcing it from a
# likelihood or inference module does not open plot devices as a side effect.
# [AI-GENERATED] Guard controlling execution when the module is sourced.
if (sys.nframe() == 0) {
  cat("\nGenerating plots...\n")
  plot_trajectory(sim)
  plot_trajectory(sim, max_time = 6, sigma_plot = 150)   # zoomed, noise visible
  plot_increment_distribution(sim)
  
  # Methods figure: the same parameter set simulated five times, showing how much
  # the trajectory varies through the randomness of the switching alone.
  plot_repeat_trajectories(v_true, lambda_true, P_true, sigma_true,
                           delta_t, T_total, n_reps = 5)
}