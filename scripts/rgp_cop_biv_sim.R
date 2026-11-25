library(copula)
library(egpd)
library(here)
library(evd)
# Simulation settings
set.seed(123)

# FOR LOCAL USE ONLY
n <- 100000       # Sample size per dataset
B <- 1000        # Number of simulation repetitions
# True data-generating copula
theta <- 2
d <- 2
cop <- gumbelCopula(theta, dim = d)

# Fixed parameters
set.seed(123)

qprob <- 0.98

xi_true <- 0.2

# Storage for results
results <- data.frame(
  iteration = seq_len(B),
  xi_hat = rep(NA_real_, B)
)

# Repeat the simulation and fitting
for (b in seq_len(B)) {

  # Generate a new dataset from the same copula
  U <- rCopula(n, cop)

  # Transform uniforms to EGPD margins
  X <- apply(
    U, 2, qegpd,
    kappa = 2,
    sigma = 1,
    xi = xi_true
  )

  R <- apply(X, 1, sum)

  quantR <- quantile(R, qprob)

  # Fit the MEGPD model
  fit <- fpot(R, quantR)

  results$xi_hat[b] <- fit$estimate[2]
  
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

# save(
#   valid,
#   summary_xi,
#   results,
#   file = here("data", paste0(scenario, ".Rdata"))
# )

# Visualize the distribution of xi estimates
# png(here("plots", paste0(scenario, ".png")))
hist(
  valid,
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

# dev.off()