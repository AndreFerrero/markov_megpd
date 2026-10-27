library(CopulaMarkov)

source("C:/Users/andre/markov_megpd/_research/r_test_scripts/megpd_2stepfit.r")

# Define the marginal distribution with its parameters
m <- egpd_margin(c(2, 1, -0.2))

# Define the dependence structure
cop <- copula::joeCopula(2)

# Construct the copula Markov model
mod <- copula_markov(m, cop)

# Simulate a time series
# (warning: Numerical inversion failures for some copulas might require multiple runs)
set.seed(46)
sim <- simulate(mod, nsim = 2000)

fit_cm <- fit(mod, init = c(2, 1, 0.2, 5), x = sim$X)

x_t <- sim$X[-1]
x_tm1 <- sim$X[-n]

R <- x_t + x_tm1
V <- log(x_t / x_tm1)

fit_results <- fit_megpd_bivariate(x_t, x_tm1, k = 10, type_s = "cr", mps = c(2, 1))

fit_radial <- fit_results$radial_fit
fit_delta <- fit_results$angular_fit

summary(fit_radial)
plot(fit_radial)

# summary(fit_delta)

# plot(fit_delta)

plot_delta(fit_delta, R)

# plot(R, V,
#      pch = 16, col = rgb(0,0,0,0.3),
#      main = "Angular spread vs radial component",
#      xlab = "R", ylab = "V")
