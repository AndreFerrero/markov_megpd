using Pkg
using MMEGPD
using Random
# using StatsPlots
using Plots

Random.seed!(1234)

n = 2000
x = simulate_markov_megpd(n, 2.0, 1.0, 0.2, :delta_exp)

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

fit_crspline = fit_markov_megpd_splines(
    x,
    K;
    spline=:crspline,
    knot_method=:quantile,
    lambda=100.0,
    x_reltol=1e-5
)

p_cr = MMEGPD.plot_spline(
    fit_crspline;
    delta_true = :delta_egpd_exp,
    show_ci = true,
    kappa = 2,
    sigma = 1,
    xi = 0.2,
    A = 0.8
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