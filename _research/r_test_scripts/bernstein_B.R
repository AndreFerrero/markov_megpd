
library(copula)
library(egpd)

set.seed(46)

# Settings
n <- 2000
nsim <- 200
m <- 8

theta <- 0.5  # Independence
xi_true <- 0.2
sigma_true <- 1
kappa_true <- 2

u_grid <- seq(0.001, 0.999, length.out = 200)

# Bernstein CDF reconstructed from weights
B_bernstein <- function(u, weights) {
  m <- length(weights)

  sapply(u, function(z) {
    if (z <= 0) return(0)
    if (z >= 1) return(1)

    sum(weights * pbeta(z, 1:m, m:1))
  })
}

# Store estimated functions
B_all <- matrix(NA, nrow = length(u_grid), ncol = nsim)

# Monte Carlo simulation
cop <- normalCopula(theta)

for (i in 1:nsim) {

  # Simulate dependent E-GPD marginals
  U <- rCopula(n, cop)

  X <- apply(
    U, 2, qegpd,
    kappa = kappa_true,
    sigma = sigma_true,
    xi = xi_true
  )

  # Aggregate
  R <- rowSums(X)

  # Fit Bernstein-EGPD
  fit <- fitegpd(
    R,
    method = "bernstein",
    bernstein.m = m
  )

  # Extract Bernstein weights
  w <- fit$bernstein.weights

  # Reconstruct B on the fitted model's scale
  B_all[, i] <- B_bernstein(u_grid, w)

  if (i %% 20 == 0) {
    cat("Replication:", i, "/", nsim, "\n")
  }
}

# Average and standard deviation
B_mean <- rowMeans(B_all, na.rm = TRUE)
B_sd <- apply(B_all, 1, sd, na.rm = TRUE)

# Plot
plot(
  u_grid, B_mean,
  type = "l",
  lwd = 2,
  ylim = c(0, 1),
  xlab = "u",
  ylab = "B(u)",
  main = paste("Average Bernstein distortion, theta =", theta)
)

# Identity reference
abline(0, 1, col = "red", lty = 2, lwd = 2)

# Monte Carlo uncertainty
lines(u_grid, pmin(1, B_mean + B_sd), lty = 3)
lines(u_grid, pmax(0, B_mean - B_sd), lty = 3)

legend(
  "topleft",
  legend = c("Average estimated B", "Identity", " +/- 1 SD"),
  col = c("black", "red", "black"),
  lty = c(1, 2, 3),
  lwd = 2,
  bty = "n"
)
