library(copula)
library(egpd)
library(here)

source(here("_research", "r_test_scripts", "megpd_2stepfit.r"))
source(here("_research","r_test_scripts","delta_copula.r"))

set.seed(46)
theta <- 4
cop <- joeCopula(theta)
xi_true <- 0.2

U <- rCopula(100000, cop)
X <- apply(U, 2, qegpd, kappa = 2, sigma = 1, xi = xi_true)
# hist(X[,2])
R <- X[, 1] + X[, 2]
# hist(R)
V <- log(X[, 1] / X[, 2])
# plot(R, V, type = "p")

fit <- fitegpd(
  R,
  method = "bernstein",
  bernstein.m = 5
)

fit$estimate

# Extract estimated weights
w <- fit$bernstein.weights
m <- length(w)

# Bernstein distortion function
B_hat <- function(u) {
  egpd:::.bernstein_cdf(u, w, m)
}

# Bernstein density
b_hat <- function(u) {
  egpd:::.bernstein_density(u, w, m)
}

# Plot
u <- seq(0, 1, length.out = 500)

plot(u, B_hat(u), type = "l", lwd = 2,
     xlab = "u", ylab = "B(u)",
     main = "Estimated Bernstein distortion")

abline(0, 1, lty = 2, col = "red")


fit <- fit_megpd_bivariate(X[, 1], X[, 2])
fit$radial_fit
plot(fit$radial_fit)

fitB <- fit_megpd_bivariate(X[, 1], X[, 2], TRUE, 10)

fitB$radial_fit
plot(fitB$radial_fit)
# fit$radial_fit


plot(hill)
# init_r <- 10^seq(-10, -1, by = 0.1)
init_r <- 0.1
r_grid <- seq(init_r, max(R), length.out = 100)

delta_est <- plot_delta(fit$angular_fit, r_grid, plot = FALSE)

# Calculate the theoretical variance at the same radius values
delta_true <- delta_copula(
  r = r_grid,
  cop = cop,
  kappa = 2,
  sigma = 1,
  xi = xi_true
)

plot(
  r_grid, delta_est$delta_hat,
  type = "l", lwd = 2, col = "black",
  main = paste0("Copula: ", getAname(cop), "; Parameter: ", theta),
  xlab = "R",
  ylab = expression(delta(R)),
  ylim = c(
    0,
    min(
      max(c(delta_est$delta_hi, delta_true)),
      3
    )
  )
)

lines(
  r_grid, delta_est$delta_lo,
  col = "red", lty = 2
)

lines(
  r_grid, delta_est$delta_hi,
  col = "red", lty = 2
)

lines(
  r_grid, delta_true,
  col = "blue", lwd = 2
)

legend(
  "topright",
  legend = c(
    "Estimated delta",
    "95% confidence interval",
    "Theoretical delta"
  ),
  col = c("black", "red", "blue"),
  lty = c(1, 2, 1),
  lwd = c(2, 1, 2),
  bty = "n"
)

rug(R, side = 1, ticksize = 0.03, col = rgb(0, 0, 0, 0.3))
