include("run_lambda_valid.jl")

using StatsBase

@load "megpd/res/lambda_eval.jld2" res_eval

function delta_curve_error(basis, beta, r_grid, delta_true_grid)

    delta_hat_fn = MMEGPD.delta_spline(basis, beta)
    delta_hat = [delta_hat_fn(r) for r in r_grid]

    return (
        rmse = sqrt(mean((delta_hat .- delta_true_grid).^2)),
        max_abs = maximum(abs.(delta_hat .- delta_true_grid)),
        curve = delta_hat
    )
end

function compute_radius(x)
    return x[1:end-1] .+ x[2:end]
end

R_train = x_train[1:end-1] .+ x_train[2:end]
R_full = data[1:end-1] .+ data[2:end]

r_grid_train = range(minimum(R_train), maximum(R_train), length=500)
r_grid_full  = range(minimum(R_full),  maximum(R_full),  length=500)

delta_true_train = [MMEGPD.get_delta(:delta_strong_upper)(r) for r in r_grid_train]
delta_true_full  = [MMEGPD.get_delta(:delta_strong_upper)(r) for r in r_grid_full]

basis_train, =
        MMEGPD.bspline_basis(
            R_train;
            nbasis=K
        )

basis_full, =
        MMEGPD.bspline_basis(
            R_full;
            nbasis=K
        )


plot_lambda = 500
plot_fit = res_eval.results[plot_lambda].fit

err_train = delta_curve_error(basis_train, plot_fit.beta, r_grid_train, delta_true_train)
err_full  = delta_curve_error(basis_full, plot_fit.beta, r_grid_full,  delta_true_full)


using Plots

p_delta = plot(
    r_grid_full,
    delta_true_full,
    linewidth = 3,
    color = :black,
    label = "true δ(r)",
    xlabel = "r",
    ylabel = "δ(r)",
    title = "δ(r): λ = $(round(plot_lambda, sigdigits=3))",
    legend = :topright
)

plot!(
    p_delta,
    r_grid_full,
    err_full.curve,
    linewidth = 2,
    color = :steelblue,
    linestyle = :solid,
    label = "estimated δ̂(r)"
)

# mark the training-support boundary — everything to the right
# of this line is extrapolation (flat, by construction of delta_spline)
vline!(
    p_delta,
    [maximum(R_train)],
    linestyle = :dash,
    color = :gray,
    linewidth = 1.5,
    label = "training support boundary"
)

display(p_delta)

println("RMSE within training support: ", err_train.rmse)
println("RMSE over full simulated range (includes extrapolation): ", err_full.rmse)

res_eval.best_lambda
best_fit = res_eval.best_fit

mean_scores = [
    res_eval.results[lambda].mean_score
    for lambda in lambda_grid
]

plot(
    lambda_grid,
    mean_scores;
    xscale = :log10,
    xlabel = "λ",
    ylabel = "Mean log predictive score",
    title = "λ selection via predictive scoring"
)

scores_distribution = [
    res_eval.results[lambda].logscores
    for lambda in lambda_grid
]

using StatsPlots

lambda_labels = string.(lambda_grid)

means = [mean(s) for s in scores_distribution]
stds  = [std(s)  for s in scores_distribution]

p = plot(
    lambda_grid,
    means;
    xscale = :log10,
    marker = :circle,
    linewidth = 2,
    xlabel = "λ",
    ylabel = "Mean log predictive score",
    title = "Predictive score versus λ",
    legend = true
)

plot!(
    p,
    lambda_grid,
    means .+ stds;
    linestyle = :dash,
    label = "± 1 SD"
)

plot!(
    p,
    lambda_grid,
    means .- stds;
    linestyle = :dash,
    label = false
)

display(p)