# MODULE 4: Diagnostics, Posterior Visualisation, and Option A vs B vs C vs D Comparison

# GENERATIVE AI ATTRIBUTION
#   Passages developed with generative AI assistance are marked inline with
#   [AI-GENERATED] or [AI-ASSISTED] tags. [AI-GENERATED] marks code
#   substantially produced by AI; [AI-ASSISTED] marks code written by me and
#   subsequently debugged, corrected or optimised with AI assistance.
#   Unmarked code is my own. Full details are in the GenAI statement.

# PURPOSE:
#   Turn the raw MCMC output from Module 3 into the figures and tables.
#   Specifically:

#     Figure 1  Trace plots (all parameters, all chains) — visual convergence.
#     Figure 2  Posterior densities with the truth overlaid, Option A vs B vs C
#               vs D on shared axes.
#     Figure 3  The sigma-recovery ladder.
#     Figure 4  Forest/caterpillar plot of posterior means +/- 95% CIs vs truth.
#     Figure 5  Running-mean (ergodic) plot for sigma — convergence over iters.
#     Table 1   Side-by-side numerical summary (A vs B vs C vs D vs truth) to CSV.

# INPUTS:
#   mcmc_A.rds, mcmc_B.rds, mcmc_C.rds, mcmc_D.rds  — saved by Module 3 (run_mcmc
#   output objects). The true parameter values, taken from Module 1's environment
#   if available, otherwise set explicitly below.

# OUTPUTS (written to ./figures/):
#   fig1_traceplots.png, fig2_posteriors.png, fig3_sigma_ladder.png,
#   fig4_forest.png, fig5_sigma_running_mean.png, and table1_summary.csv.


# SECTION 0: Load cached runs and set up

# Load the two posterior objects. If they are already in the workspace (e.g.
# already present in the workspace), those are used; otherwise they are read
# from disk.
if (!exists("mcmc_A")) mcmc_A <- readRDS("mcmc_A.rds")
if (!exists("mcmc_B")) mcmc_B <- readRDS("mcmc_B.rds")
if (!exists("mcmc_C")) mcmc_C <- readRDS("mcmc_C.rds")
if (!exists("mcmc_D")) mcmc_D <- readRDS("mcmc_D.rds")

# True parameter values for overlay. Prefer the live values from Module 1;
# fall back to the known simulation truth if Module 1 is not sourced.
if (!exists("v_true"))      v_true      <- c(2000, -1500)
if (!exists("lambda_true")) lambda_true <- c(1.0, 0.5)
if (!exists("sigma_true"))  sigma_true  <- 50

# Truth on the sampling (theta) scale: c(v1, v2, log_lambda1, log_lambda2, log_sigma)
theta_true <- c(v_true, log(lambda_true), log(sigma_true))
# Truth on the natural scale: c(v1, v2, lambda1, lambda2, sigma)
nat_true   <- c(v_true, lambda_true, sigma_true)

param_names    <- mcmc_A$param_names                      # theta-scale names
nat_names      <- c(paste0("v", seq_along(v_true)),
                    paste0("lambda", seq_along(lambda_true)), "sigma")
n_states       <- mcmc_A$prior$n_states
d              <- length(param_names)

# Make sure an output directory exists
dir.create("figures", showWarnings = FALSE)

# Pull the post-burn-in samples of one chain as a matrix
get_post <- function(mcmc, chain_idx) {
  post_idx <- (mcmc$burn_in + 1):mcmc$n_iter
  mcmc$chains[[chain_idx]]$chain[post_idx, , drop = FALSE]
}

# Self-contained Gelman-Rubin R-hat (a copy of Module 3's, so Module 4 can run
# from the cached .rds files without re-sourcing Module 3 and triggering the
# ~15-minute sampler). Identical formula to Module 3.
gelman_rubin <- function(mcmc_out) {
  chains  <- mcmc_out$chains
  burn_in <- mcmc_out$burn_in
  n_iter  <- mcmc_out$n_iter
  m       <- mcmc_out$n_chains
  d       <- length(mcmc_out$param_names)
  post_idx <- (burn_in + 1):n_iter
  n        <- length(post_idx)
  rhat <- numeric(d)
  for (p in 1:d) {
    M <- sapply(chains, function(ch) ch$chain[post_idx, p])
    chain_means <- colMeans(M); grand_mean <- mean(chain_means)
    B <- (n / (m - 1)) * sum((chain_means - grand_mean)^2)
    W <- mean(apply(M, 2, var))
    var_hat <- ((n - 1) / n) * W + (1 / n) * B
    rhat[p] <- sqrt(var_hat / W)
  }
  names(rhat) <- mcmc_out$param_names
  rhat
}

# Pooled post-burn-in samples (all chains stacked), theta scale
pool_theta <- function(mcmc) {
  do.call(rbind, lapply(seq_len(mcmc$n_chains), function(c) get_post(mcmc, c)))
}

# Convert a theta-scale sample matrix to the natural scale (exp the last n+1 cols)
to_natural <- function(M, n) {
  out <- M
  cols <- c((n + 1):(2 * n), 2 * n + 1)   # log_lambda... and log_sigma
  out[, cols] <- exp(M[, cols])
  out
}

col_A   <- "orange2"      # Option A
col_B   <- "dodgerblue3"  # Option B
col_C   <- "seagreen4"    # Option C
col_D   <- "firebrick3"   # Option D
col_tru <- "black"        # truth

#   A = solid (1), B = dashed (2), C = dotted (3), D = dot-dash (4).
lty_A <- 1; lty_B <- 2; lty_C <- 3; lty_D <- 4

# [AI-ASSISTED] Chain colour palette. An earlier version indexed the palette
# such that chains 5 and 6 of a six-chain run were drawn in a fully
# transparent colour and were invisible; this was identified and corrected
# with AI assistance. The palette below was chosen with AI assistance to be
# distinguishable with colour-blindness, and is paired with distinct line
# types so that the figures remain readable in black and white.
# See GenAI statement.
col_chain <- c("dodgerblue3", "orange2", "seagreen4", "orchid3",
               "firebrick3", "goldenrod3")  # up to 6 chains


# SECTION 1: Figure 1 — Trace plots (convergence check)

# One panel per parameter; the four chains overlaid in different colours.
# Well-mixed, converged chains look like overlapping "fuzzy caterpillars"
# with no trends or separated bands. We show Option B; change `mc` to mcmc_A to
# produce the Option-A version.

make_traceplots <- function(mc, tag) {
  png(sprintf("figures/fig1_traceplots_%s.png", tag),
      width = 1100, height = 1400, res = 130)
  op <- par(mfrow = c(d, 1), mar = c(2.4, 4, 1.6, 1), oma = c(2, 0, 2, 0))
  
  post_idx <- (mc$burn_in + 1):mc$n_iter
  for (p in 1:d) {
    # y-range across all chains for this parameter
    yr <- range(sapply(seq_len(mc$n_chains),
                       function(c) range(mc$chains[[c]]$chain[post_idx, p])))
    plot(NA, xlim = c(1, length(post_idx)), ylim = yr,
         xlab = "", ylab = param_names[p], cex.lab = 1.1)
    for (c in seq_len(mc$n_chains)) {
      lines(mc$chains[[c]]$chain[post_idx, p],
            col = adjustcolor(col_chain[c], alpha.f = 0.7), lwd = 0.5)
    }
    # truth line (theta scale)
    abline(h = theta_true[p], col = col_tru, lwd = 2, lty = 2)
  }
  mtext(sprintf("Trace plots — Option %s (%d chains, post burn-in)",
                tag, mc$n_chains),
        outer = TRUE, cex = 1.1, font = 2)
  mtext("iteration (post burn-in)", side = 1, outer = TRUE)
  par(op); dev.off()
}


# SECTION 2: Figure 2 — Posterior densities, Option A vs B vs C vs D, truth overlaid

# For each parameter on the natural scale, plot the kernel-density estimate of
# the pooled posterior under Option A, Option B and Option C on shared axes,
# with a vertical line at the true value. This is the single most informative
# figure: Option C and Option B's densities should sit closer to the truth
# (especially sigma and the velocities) than Option A's.

make_posteriors <- function() {
  natA <- to_natural(pool_theta(mcmc_A), n_states); colnames(natA) <- nat_names
  natB <- to_natural(pool_theta(mcmc_B), n_states); colnames(natB) <- nat_names
  natC <- to_natural(pool_theta(mcmc_C), n_states); colnames(natC) <- nat_names
  natD <- to_natural(pool_theta(mcmc_D), n_states); colnames(natD) <- nat_names
  
  png("figures/fig2_posteriors.png", width = 1300, height = 850, res = 130)
  op <- par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1))
  
  for (p in 1:d) {
    dA <- density(natA[, p]); dB <- density(natB[, p])
    dC <- density(natC[, p]); dD <- density(natD[, p])
    xr <- range(dA$x, dB$x, dC$x, dD$x, nat_true[p])
    yr <- c(0, max(dA$y, dB$y, dC$y, dD$y) * 1.05)
    plot(NA, xlim = xr, ylim = yr, main = nat_names[p],
         xlab = nat_names[p], ylab = "posterior density")
    polygon(dA$x, dA$y, col = adjustcolor(col_A, 0.22), border = col_A,
            lwd = 2, lty = lty_A)
    polygon(dB$x, dB$y, col = adjustcolor(col_B, 0.22), border = col_B,
            lwd = 2, lty = lty_B)
    polygon(dC$x, dC$y, col = adjustcolor(col_C, 0.22), border = col_C,
            lwd = 2, lty = lty_C)
    polygon(dD$x, dD$y, col = adjustcolor(col_D, 0.22), border = col_D,
            lwd = 2, lty = lty_D)
    abline(v = nat_true[p], col = col_tru, lwd = 2, lty = 2)
  }
  # shared legend in the 6th panel space
  plot.new()
  legend("center",
         legend = c("Option A (zero-switch)", "Option B (one-switch)",
                    "Option C (two-switch)", "Option D (exact)", "truth"),
         col = c(col_A, col_B, col_C, col_D, col_tru), lwd = 2,
         lty = c(lty_A, lty_B, lty_C, lty_D, 2),
         bty = "n", cex = 1.1)
  par(op); dev.off()
}


# SECTION 3: Figure 3 — sigma-recovery across the model refinements

# Shows how the posterior estimate of sigma moves toward the true value of 50
# as the model is refined. The first ("sigma^2 bug") value is the pre-fix
# Option-A result recorded during development; the rest come from the cached
# runs.

make_sigma_ladder <- function() {
  # Pull sigma posteriors (natural scale)
  sigA <- exp(pool_theta(mcmc_A)[, 2 * n_states + 1])
  sigB <- exp(pool_theta(mcmc_B)[, 2 * n_states + 1])
  sigC <- exp(pool_theta(mcmc_C)[, 2 * n_states + 1])
  sigD <- exp(pool_theta(mcmc_D)[, 2 * n_states + 1])
  
  stages <- c("Option A\n(sigma^2 bug)*",
              "Option A\n(2*sigma^2 fix)",
              "Option B\n(one-switch)",
              "Option C\n(two-switch)",
              "Option D\n(exact)",
              "Truth")
  # NOTE ON PROVENANCE OF THE FIRST RUNG (148.3):
  # Every other value on this ladder is computed live from the cached MCMC
  # runs (sigA/sigB/sigC, from mcmc_A/B/C.rds). The first rung, 148.3, is the
  # only hard-coded number. It is the Option-A posterior mean for sigma from an
  # early development run in which the measurement-noise variance was wrongly
  # set to sigma^2 instead of 2*sigma^2 (the variance of the differenced
  # observation noise). After that fix, Option A gives ~105 (the second rung),
  # so 148.3 cannot be reproduced by the current pipeline; it is retained only
  # to document the full diagnostic progression of the project. It has no
  # credible interval (los/his = NA below) for the same reason.
  means  <- c(148.3, mean(sigA), mean(sigB), mean(sigC), mean(sigD), sigma_true)
  los    <- c(NA, quantile(sigA, .025), quantile(sigB, .025),
              quantile(sigC, .025), quantile(sigD, .025), NA)
  his    <- c(NA, quantile(sigA, .975), quantile(sigB, .975),
              quantile(sigC, .975), quantile(sigD, .975), NA)
  cols   <- c(col_A, col_A, col_B, col_C, col_D, col_tru)
  # Distinct plot symbols per option so the rungs are also distinguishable in
  # black and white: both Option-A rungs = circle (19), B = triangle (17),
  # C = square (15), D = diamond (18), Truth = plus (3).
  pchs   <- c(19, 19, 17, 15, 18, 3)
  
  png("figures/fig3_sigma_ladder.png", width = 1100, height = 720, res = 130)
  op <- par(mar = c(5, 5, 3, 1))
  x <- seq_along(stages)
  plot(x, means, ylim = range(40, 160), pch = pchs, cex = 1.8, col = cols,
       xaxt = "n", xlab = "", ylab = expression(sigma~"(posterior mean, 95% CI)"),
       main = expression("Recovery of the noise parameter "*sigma*" across model refinements"))
  axis(1, at = x, labels = stages, cex.axis = 0.8, padj = 0.5)
  arrows(x, los, x, his, angle = 90, code = 3, length = 0.06, col = cols, lwd = 2)
  abline(h = sigma_true, col = col_tru, lwd = 1.5, lty = 2)
  text(x, means, labels = sprintf("%.1f", means), pos = 3, cex = 0.85, offset = 0.8)
  mtext("* pre-fix value recorded during development (not from cached run)",
        side = 1, line = 3.4, cex = 0.7, adj = 0)
  par(op); dev.off()
}


# SECTION 4: Figure 4 — Forest plot (posterior mean +/- 95% CI vs truth)

# A compact "did we recover the truth?" figure: for every parameter, a point
# at the posterior mean with a 95% CI whisker, Option A above Option B, with
# the true value marked. Parameters are shown on the natural scale but each in
# its own panel because their scales differ by orders of magnitude.

make_forest <- function() {
  natA <- to_natural(pool_theta(mcmc_A), n_states); colnames(natA) <- nat_names
  natB <- to_natural(pool_theta(mcmc_B), n_states); colnames(natB) <- nat_names
  natC <- to_natural(pool_theta(mcmc_C), n_states); colnames(natC) <- nat_names
  natD <- to_natural(pool_theta(mcmc_D), n_states); colnames(natD) <- nat_names
  
  png("figures/fig4_forest.png", width = 1300, height = 850, res = 130)
  op <- par(mfrow = c(2, 3), mar = c(4, 2, 2.5, 1))
  
  for (p in 1:d) {
    mA <- mean(natA[, p]); ciA <- quantile(natA[, p], c(.025, .975))
    mB <- mean(natB[, p]); ciB <- quantile(natB[, p], c(.025, .975))
    mC <- mean(natC[, p]); ciC <- quantile(natC[, p], c(.025, .975))
    mD <- mean(natD[, p]); ciD <- quantile(natD[, p], c(.025, .975))
    xr <- range(ciA, ciB, ciC, ciD, nat_true[p])
    plot(NA, xlim = xr, ylim = c(0.5, 4.5), yaxt = "n",
         main = nat_names[p], xlab = nat_names[p], ylab = "")
    axis(2, at = c(4, 3, 2, 1), labels = c("A", "B", "C", "D"), las = 1)
    abline(v = nat_true[p], col = col_tru, lwd = 2, lty = 2)
    # Option A at y=4 (circle)
    arrows(ciA[1], 4, ciA[2], 4, angle = 90, code = 3, length = 0.05,
           col = col_A, lwd = 2, lty = lty_A)
    points(mA, 4, pch = 19, col = col_A, cex = 1.4)
    # Option B at y=3 (triangle)
    arrows(ciB[1], 3, ciB[2], 3, angle = 90, code = 3, length = 0.05,
           col = col_B, lwd = 2, lty = lty_B)
    points(mB, 3, pch = 17, col = col_B, cex = 1.4)
    # Option C at y=2 (square)
    arrows(ciC[1], 2, ciC[2], 2, angle = 90, code = 3, length = 0.05,
           col = col_C, lwd = 2, lty = lty_C)
    points(mC, 2, pch = 15, col = col_C, cex = 1.4)
    # Option D at y=1 (diamond)
    arrows(ciD[1], 1, ciD[2], 1, angle = 90, code = 3, length = 0.05,
           col = col_D, lwd = 2, lty = lty_D)
    points(mD, 1, pch = 18, col = col_D, cex = 1.6)
  }
  plot.new()
  legend("center",
         legend = c("Option A", "Option B", "Option C", "Option D", "truth"),
         col = c(col_A, col_B, col_C, col_D, col_tru), lwd = 2,
         lty = c(lty_A, lty_B, lty_C, lty_D, 2),
         pch = c(19, 17, 15, 18, NA), bty = "n", cex = 1.1)
  par(op); dev.off()
}


# SECTION 5: Figure 5 — Running mean of sigma (ergodic convergence)

# The running (cumulative) mean of sigma for each chain across post-burn-in
# iterations. Flat, overlapping lines converging to a common value are visual
# evidence the chains have equilibrated. Shown for Option B.

make_running_mean <- function(mc, tag) {
  sigma_col <- 2 * n_states + 1
  post_idx  <- (mc$burn_in + 1):mc$n_iter
  
  png(sprintf("figures/fig5_sigma_running_mean_%s.png", tag),
      width = 1000, height = 640, res = 130)
  op <- par(mar = c(4.5, 4.5, 2.5, 1))
  plot(NA, xlim = c(1, length(post_idx)), ylim = c(45, 75),
       xlab = "iteration (post burn-in)",
       ylab = expression("running mean of "*sigma),
       main = sprintf("Ergodic mean of sigma — Option %s", tag))
  for (c in seq_len(mc$n_chains)) {
    s <- exp(mc$chains[[c]]$chain[post_idx, sigma_col])
    lines(cumsum(s) / seq_along(s), col = col_chain[c], lwd = 1.5)
  }
  abline(h = sigma_true, col = col_tru, lwd = 2, lty = 2)
  legend("topright", legend = paste("chain", seq_len(mc$n_chains)),
         col = col_chain[seq_len(mc$n_chains)], lwd = 1.5, bty = "n", cex = 0.9)
  par(op); dev.off()
}


# SECTION 6: Table 1 — Numerical summary (A vs B vs C vs D vs truth) to CSV

make_summary_table <- function() {
  natA <- to_natural(pool_theta(mcmc_A), n_states); colnames(natA) <- nat_names
  natB <- to_natural(pool_theta(mcmc_B), n_states); colnames(natB) <- nat_names
  natC <- to_natural(pool_theta(mcmc_C), n_states); colnames(natC) <- nat_names
  natD <- to_natural(pool_theta(mcmc_D), n_states); colnames(natD) <- nat_names
  rhatA <- gelman_rubin(mcmc_A)
  rhatB <- gelman_rubin(mcmc_B)
  rhatC <- gelman_rubin(mcmc_C)
  rhatD <- gelman_rubin(mcmc_D)
  
  tab <- data.frame(
    parameter = nat_names,
    truth     = nat_true,
    A_mean    = apply(natA, 2, mean),
    A_q2.5    = apply(natA, 2, quantile, .025),
    A_q97.5   = apply(natA, 2, quantile, .975),
    A_rhat    = as.numeric(rhatA),
    B_mean    = apply(natB, 2, mean),
    B_q2.5    = apply(natB, 2, quantile, .025),
    B_q97.5   = apply(natB, 2, quantile, .975),
    B_rhat    = as.numeric(rhatB),
    C_mean    = apply(natC, 2, mean),
    C_q2.5    = apply(natC, 2, quantile, .025),
    C_q97.5   = apply(natC, 2, quantile, .975),
    C_rhat    = as.numeric(rhatC),
    D_mean    = apply(natD, 2, mean),
    D_q2.5    = apply(natD, 2, quantile, .025),
    D_q97.5   = apply(natD, 2, quantile, .975),
    D_rhat    = as.numeric(rhatD),
    row.names = NULL
  )
  # percentage error of the posterior mean relative to truth
  tab$A_pct_err <- round(100 * (tab$A_mean - tab$truth) / tab$truth, 1)
  tab$B_pct_err <- round(100 * (tab$B_mean - tab$truth) / tab$truth, 1)
  tab$C_pct_err <- round(100 * (tab$C_mean - tab$truth) / tab$truth, 1)
  tab$D_pct_err <- round(100 * (tab$D_mean - tab$truth) / tab$truth, 1)
  
  write.csv(tab, "figures/table1_summary.csv", row.names = FALSE)
  tab
}


# SECTION 7: Build everything

cat("\n=== MODULE 4: DIAGNOSTICS & COMPARISON FIGURES ===\n")

make_traceplots(mcmc_D, "D")
make_traceplots(mcmc_C, "C")
make_traceplots(mcmc_B, "B")
make_traceplots(mcmc_A, "A")
cat("  [1/6] trace plots written (Options A, B, C, D)\n")

make_posteriors()
cat("  [2/6] posterior-density comparison written\n")

make_sigma_ladder()
cat("  [3/6] sigma-recovery ladder written\n")

make_forest()
cat("  [4/6] forest plot written\n")

make_running_mean(mcmc_B, "B")
cat("  [5/6] sigma running-mean plot written\n")

summary_table <- make_summary_table()
cat("  [6/6] summary table written\n\n")

cat("Summary table (natural scale):\n")
print(summary_table, digits = 4, row.names = FALSE)

cat("\nAll figures saved to ./figures/.\n")