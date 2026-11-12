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

library(JuliaCall)
julia_setup()

julia_source("megpd/code/src/MMEGPD.jl")

K <- as.integer(5)
fit <- julia_call(
    "MMEGPD.fit_markov_megpd_splines",
    egpd_gumbel_data$x,
    K
)

fit$result
fit$theta

r_grid <- seq(0, max(R), length.out = 500)
julia_eval("using BSplines")
basis <- julia_call("MMEGPD.bspline_basis", R, nbasis = K, degree = as.integer(3))
B_grid <- julia_call("basismatrix", basis[[1]], r_grid)

log_delta_hat <- B_grid %*% fit$beta
delta_hat <- exp(log_delta_hat)

plot(r_grid, delta_hat, type = "l")
