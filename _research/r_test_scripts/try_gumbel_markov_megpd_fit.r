source("megpd/code/megpd_fit.r")
source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/copulas/gumbel.r")

egpd_gumbel_model <- make_copula_markov_model(margin_egpd, copula_gumbel, stan_mod = NULL)

set.seed(123)
margin_param <- c(mu = 0, kappa = 2, sigma = 1, xi = 0.2)
copula_param <- 4
n <- 6000

egpd_gumbel_data <- egpd_gumbel_model$simulate(
    n = n,
    margin_param = margin_param,
    copula_param = copula_param
)

x_t <- egpd_gumbel_data$x[-1]
x_tm1 <- egpd_gumbel_data$x[-n]

R <- x_t + x_tm1
V <- log(x_t / x_tm1)

fit_results <- fit_megpd_bivariate(x_t, x_tm1, k = 10, type_s = "cr", mps = c(2, 1))

fit_radial <- fit_results$radial_fit
fit_delta <- fit_results$angular_fit

# summary(fit_radial)
# plot(fit_radial)

# summary(fit_delta)

# plot(fit_delta)

plot_delta(fit_delta, R)

# plot(R, V,
#      pch = 16, col = rgb(0,0,0,0.3),
#      main = "Angular spread vs radial component",
#      xlab = "R", ylab = "V")
