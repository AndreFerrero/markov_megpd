using DrWatson
@quickactivate "markov_megpd"

using Pkg
using MMEGPD
using Random
# using StatsPlots
using Plots

Random.seed!(1234)

n = 5000
theta0 = [2.0, 1.0, 0.2]
x = MMEGPD.simulate_markov_megpd(n, theta0, :delta_exp; x0 = 30)

R = x[1:(end-1)] .+ x[2:end]
V = x[2:end] ./ x[1:end-1]

r_grid = range(
    minimum(R),
    maximum(R),
    length=500
)

K = 10

# fit_pspline = MMEGPD.fit_markov_megpd_splines(
#     x,
#     K;
#     spline=:pspline,
#     lambda=8.73,
#     p=2,
#     x_reltol=1e-6
# )

# p_ps = MMEGPD.plot_spline(
#     fit_pspline;
#     delta_true = :delta_strong_upper
# )

# # CR-spline knots
# vline!(
#     p_ps,
#     fit_pspline.basis.breakpoints;
#     label = "B-knots",
#     alpha = 0.5,
#     linestyle = :dash
# )

# # Original R locations used to construct the spline
# scatter!(
#     p_ps,
#     R,
#     fill(first(ylims(p_ps)), length(R));
#     marker = :vline,
#     markersize = 7,
#     label = "R points"
# )

# display(p_ps)

fit_crspline = MMEGPD.fit_markov_megpd_splines(
    x,
    K;
    bernstein=false,
    lambda=100.0,
    x_reltol=1e-4,
    verbose=true
)

fit_crspline.theta
fit_crspline.w
fit_crspline.beta

MMEGPD.radial_qqplot(fit_crspline)

p_cr = MMEGPD.plot_spline(
    fit_crspline;
    delta_true = :delta_exp,
    show_ci = true
)

# CR-spline knots
vline!(
    p_cr,
    fit_crspline.basis.breakpoints;
    label = "CR-knots",
    alpha = 0.5,
    linestyle = :dash
)

# Original R locations used to construct the spline
scatter!(
    p_cr,
    fit_crspline.R,
    fill(first(ylims(p_cr)), length(fit_crspline.R));
    marker = :vline,
    markersize = 7,
    label = "R points"
)

display(p_cr)