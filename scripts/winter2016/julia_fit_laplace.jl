
using MMEGPD
using RCall
using StatsPlots
using StatsBase

R"source(\"winter2016/load_data.R\")"

x = rcopy(R"j.data")

R = x[:,1] .+ x[:,2]
V = x[:,1] ./ x[:,2]

r_grid =
    range(
        minimum(R),
        maximum(R),
        length=500
    )

scatter(
    R,
    log.(V);
    xlabel = "R",
    ylabel = "log(V)",
    label = "Data",
    alpha = 0.4,
    legend = :topright,
)

histogram(R)

K = 10

test_lambda = 10
test_p = 2

test_fit = MMEGPD.fit_markov_megpd_splines(
    x, K;
    degree=3,
    p=test_p,
    lambda=test_lambda,
    map=false,
    verbose=true,
    x_reltol=1e-5,
    kappa_init=30,
    sigma_init=10
)

test_fit.theta

MMEGPD.plot_spline(test_fit, r_grid; plot_fn = true)
plot!(title = "λ = $test_lambda; p = $test_p", legend = :topright)
plot!(ylims = (0, 0.3))


grid_plot = MMEGPD.evaluate_lambda_grid(
    x, K; p = 2, lower = 0.0, ngrid = 15, x_reltol = 1e-5, verbose = true
)

using Plots
plot(
    grid_plot.log10lambda,
    grid_plot.logL,
    marker=:circle,
    xlabel="log10(λ)",
    ylabel="Laplace marginal log likelihood",
    title="Laplace objective; K = $K, p = 2",
    legend=false
)

plot(
    grid_plot.lambda,
    grid_plot.logL,
    marker=:circle,
    xlabel="λ",
    ylabel="Laplace marginal log likelihood",
    title="Laplace objective; K = $K, p = 2",
    legend=false
)

i = argmax(grid_plot.logL)

max_logL = grid_plot.logL[i]
max_lambda = grid_plot.lambda[i]
max_log10lambda = grid_plot.log10lambda[i]

println("Maximum logL = ", max_logL)
println("λ = ", max_lambda)
println("log10(λ) = ", max_log10lambda)

la_res_optimonly = MMEGPD.laplace_optimize(x,K; p = 2, map = true, verbose = true, x_reltol=1e-5)
