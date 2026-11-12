using DrWatson
@quickactivate "markov_megpd"

using RCall
using MMEGPD
using Plots

R"""
library(CopulaMarkov)

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
"""

x = rcopy(R"sim$X")

K = 10
# fit = MMEGPD.fit_markov_megpd_splines(
#     x, K;
#     x_reltol = 1e-6,
#     lambda = 50,
#     sigma_init = 3,
#     xi_init = -0.2
# )

# fit.theta

la_fit = MMEGPD.laplace_optimize(
    x,
    K;
    spline=:crspline,
    knot_method=:quantile,
    x_reltol=1e-5,
    verbose=true,
    sigma_init = 3.0,
    xi_init = -0.2
)

fit = la_fit.best_fit

fit.theta
fit.lambda

p = MMEGPD.plot_spline(
    fit
)
title!("λ = $(round(fit.lambda; digits = 2))")
# CR-spline knots
vline!(
    p,
    fit.basis.breakpoints;
    label = "knots",
    alpha = 0.5,
    linestyle = :dash
)

# Original R locations used to construct the spline
scatter!(
    p,
    fit.R,
    fill(first(ylims(p)), length(fit.R));
    marker = :vline,
    markersize = 7,
    label = "R points"
)

display(p)
