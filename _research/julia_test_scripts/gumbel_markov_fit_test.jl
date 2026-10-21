using RCall

using MMEGPD

x = MMEGPD.sample_gumbel_markov(2000, 2, 1, 0.2, 3)

R = x[1:(end-1)] .+ x[2:end]
V = x[2:end] ./ x[1:end-1]

r_grid = range(minimum(R), maximum(R), 500)

K = 8

fit_pspline = MMEGPD.fit_markov_megpd_splines(
    x,
    K;
    spline=:pspline,
    lambda=10,
    p=2,
    x_reltol=1e-6
)

fit_crspline = MMEGPD.fit_markov_megpd_splines(
    x,
    K;
    spline=:crspline,
    knot_method=:quantile,
    lambda=1145,
    x_reltol=1e-6
)


p_ps = MMEGPD.plot_spline(
    fit_pspline
)

# CR-spline knots
vline!(
    p_ps,
    fit_pspline.basis.breakpoints;
    label = "B-knots",
    alpha = 0.5,
    linestyle = :dash
)

# Original R locations used to construct the spline
scatter!(
    p_ps,
    R,
    fill(first(ylims(p_ps)), length(R));
    marker = :vline,
    markersize = 7,
    label = "R points"
)

display(p_ps)

p_cr = MMEGPD.plot_spline(
    fit_crspline, r_grid
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
    R,
    fill(first(ylims(p_cr)), length(R));
    marker = :vline,
    markersize = 7,
    label = "R points"
)

display(p_cr)