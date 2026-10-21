using MMEGPD
using Random
using Plots

Random.seed!(1234)


x = MMEGPD.simulate_markov_megpd(2000, 2, 1, 0.2, :delta_strong_upper)
# x = MMEGPD.sample_gumbel_markov(2000, 2, 1, 0.2, 3)


R = x[1:(end-1)] .+ x[2:end]

K = 10

# Log marginal likelihood evaluation

# grid_plot = MMEGPD.laplace_grid_search(
#     x, K;
#     spline=:crspline, knot_method=:even,
#     ngrid=18,
#     lower=0.0, upper=4.0,
#     x_reltol=1e-6, verbose=true
# )

# plot(
#     grid_plot.log10lambda,
#     grid_plot.logL,
#     marker=:circle,
#     xlabel="log10(λ)",
#     ylabel="Laplace marginal log likelihood",
#     title="Laplace objective; K = $K",
#     legend=false
# )

# plot(
#     grid_plot.lambda,
#     grid_plot.logL,
#     marker=:circle,
#     xlabel="λ",
#     ylabel="Laplace marginal log likelihood",
#     title="Laplace objective; K = $K",
#     legend=false
# )

# i = argmax(grid_plot.logL)

# max_logL = grid_plot.logL[i]
# max_lambda = grid_plot.lambda[i]
# max_log10lambda = grid_plot.log10lambda[i]

# println("Maximum logL = ", max_logL)
# println("λ = ", max_lambda)
# println("log10(λ) = ", max_log10lambda)

# Optimisation routine for lambda selection
la_res_optimonly = MMEGPD.laplace_optimize(
    x, K;
    spline=:crspline, knot_method=:quantile,
    map=true, hyper_inv_sin = true,
    lower=-2.0, upper=2.0,
    verbose=true, x_reltol=1e-6
)

# Visualisation of delta function

plot_fit = la_res_optimonly.best_fit
plot_lambda = la_res_optimonly.best_lambda

plot_fit.theta

plot_delta = MMEGPD.plot_spline(plot_fit;
    delta_true = :delta_strong_upper,
    xscale = :asinh
)
plot!(title="λ = $(round(plot_lambda, sigdigits=3))")
# plot!(ylims=(0, 1))

# CR-spline knots
vline!(
    plot_fit.basis.breakpoints;
    label = "CR-knots",
    alpha = 0.5,
    linestyle = :dash
)

# Original R locations used to construct the spline
scatter!(
    plot_fit.R,
    fill(first(ylims(plot_delta)), length(R));
    marker = :vline,
    markersize = 7,
    label = "R points"
)