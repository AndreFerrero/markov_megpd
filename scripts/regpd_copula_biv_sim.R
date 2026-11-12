library(copula)
library(egpd)
library(here)

# Simulation settings
set.seed(123)

n <- 100000       # Sample size per dataset
B <- 100        # Number of simulation repetitions
xi_true <- 0.2
# True data-generating copula
theta <- 2
cop <- joeCopula(theta)
m_degree <- 10

# Storage for results
results <- data.frame(
  iteration = seq_len(B),
  xi_hat = rep(NA_real_, B),
  sigma_hat = rep(NA_real_, B),
  kappa_hat = rep(NA_real_, B),
  converged = rep(FALSE, B)
)

# Repeat the simulation and fitting
for (b in seq_len(B)) {

  # Generate a new dataset from the same Joe copula
  U <- rCopula(n, cop)

  # Transform uniforms to EGPD margins
  X <- apply(
    U, 2, qegpd,
    kappa = 2,
    sigma = 1,
    xi = xi_true
  )

  R <- X[, 1] + X[, 2]
  # Fit the MEGPD model
  fit <- tryCatch(
    fitegpd(
      R,
      method = "bernstein",
      bernstein.m = m_degree
      ),
    error = function(e) {
      message("Iteration ", b, " failed: ", conditionMessage(e))
      NULL
    }
  )

  # Store estimates if fitting succeeded
  if (!is.null(fit)) {

    radial_est <- fit$estimate

    if (!is.null(radial_est) && "xi" %in% names(radial_est)) {
      results$xi_hat[b] <- unname(radial_est["xi"])
      results$converged[b] <- is.finite(results$xi_hat[b])
    }

    if (!is.null(radial_est) && "sigma" %in% names(radial_est)) {
      results$sigma_hat[b] <- unname(radial_est["sigma"])
    }

    if (!is.null(radial_est) && "kappa" %in% names(radial_est)) {
      results$kappa_hat[b] <- unname(radial_est["kappa"])
    }
  }

  if (b %% 5 == 0) {
    message("Completed ", b, " / ", B)
  }
}

# Keep successful estimates
valid <- results$xi_hat[is.finite(results$xi_hat)]

# Summarize xi estimates
summary_xi <- c(
  repetitions = B,
  successful_fits = length(valid),
  mean_xi = mean(valid),
  bias_vs_margin_xi = mean(valid) - (xi_true),
  sd_xi = sd(valid),
  median_xi = median(valid),
  q025 = unname(quantile(valid, 0.025)),
  q975 = unname(quantile(valid, 0.975)),
  min_xi = min(valid),
  max_xi = max(valid)
)

print(summary_xi)

# Inspect individual estimates
print(head(results))

save(valid, summary_xi,
file = here("data", "joe2_radialegpd_m10B100n100k.Rdata"))

# Visualize the distribution of xi estimates
png(here("plots", "joe2_radialegpd_m10B100n100k.png"))
hist(
  valid,
  breaks = 30,
  xlab = expression(hat(xi)),
  main = paste0("Copula: ", getAname(cop), "; Parameter: ", theta),
  col = "lightblue",
  border = "white"
)

abline(v = mean(valid), col = "red", lwd = 2)
abline(v = xi_true, col = "darkgreen", lwd = 2, lty = 2)

legend(
  "topright",
  legend = c("Mean estimate", paste0("Marginal xi = ", xi_true)),
  col = c("red", "darkgreen"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

dev.off()