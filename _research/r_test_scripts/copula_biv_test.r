library(copula)
library(egpd)
library(here)

source(here("_research", "r_test_scripts", "megpd_2stepfit.r"))
source(here("_research","r_test_scripts","delta_copula.r"))

set.seed(46)
theta <- 6
cop <- joeCopula(theta)
xi_true <- -0.2

U <- rCopula(2000, cop)
X <- apply(U, 2, qegpd, kappa = 2, sigma = 1, xi = xi_true)
# hist(X[,2])
R <- X[, 1] + X[, 2]
# hist(R)
V <- log(X[, 1] / X[, 2])
# plot(R, V, type = "p")

fit <- fit_megpd_bivariate(X[, 1], X[, 2])
# fit$radial_fit

# init_r <- 10^seq(-10, -1, by = 0.1)
init_r <- 0.01
after_r <- seq(init_r[length(init_r)], max(R), length.out = 100)
r_grid <- c(init_r, after_r)

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
