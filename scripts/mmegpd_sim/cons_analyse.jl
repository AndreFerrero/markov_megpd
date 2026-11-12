using DrWatson
@quickactivate "markov_megpd"

save_file = joinpath(datadir(), "cons_bern_egpd_crs_nrep100_n2k_K5_xi02_deltaexp.jld2")

using JLD2
using Statistics
using DataFrames
using StatsPlots
using MMEGPD
using Optim

@load save_file fits nrep n K kappa xi delta_true

# Extract theta = (kappa, xi) from each fitted model
theta_estimates = [fit.theta for fit in fits]

# Extract only kappa and xi
theta_matrix = reduce(vcat, [
    collect(theta)' for theta in theta_estimates
])

# Create DataFrame
df = DataFrame(
    replicate = 1:length(fits),
    kappa_hat = theta_matrix[:, 1],
    xi_hat = theta_matrix[:, 2]
)

println(df)

mean(df[:,2])
median(df[:,2])

p1 = @df df boxplot(
    :kappa_hat,
    title = "Estimation of κ",
    label = "Estimates",
    color = :skyblue,
    alpha = 0.7,
    linewidth = 1.5
)

hline!(
    p1,
    [kappa],
    color = :red,
    linestyle = :dash,
    linewidth = 2,
    label = "True value"
)

p2 = @df df boxplot(
    :xi_hat,
    title = "Estimation of ξ",
    label = "Estimates",
    color = :lightgreen,
    alpha = 0.7,
    linewidth = 1.5
)

hline!(
    p2,
    [xi],
    color = :red,
    linestyle = :dash,
    linewidth = 2,
    label = "True value"
)

p_box = plot(
    p1, p2,
    layout = (1, 2),
    size = (1000, 450)
)

display(p_box)



MMEGPD.plot_spline(fits[30]; delta_true=:delta_exp)


# Average Bernstein cdf

# Grid
u_grid = range(0, 1, length=200)

# Evaluate each Bernstein CDF
B_values = [
    [MMEGPD.bernstein_cdf(u, fit.w) for u in u_grid]
    for fit in fits
]

# Average CDF
B_mean = [
    mean(B_values[j][i] for j in eachindex(fits))
    for i in eachindex(u_grid)
]

# Pointwise 95% confidence bands
B_lower = [
    quantile([B_values[j][i] for j in eachindex(fits)], 0.025)
    for i in eachindex(u_grid)
]

B_upper = [
    quantile([B_values[j][i] for j in eachindex(fits)], 0.975)
    for i in eachindex(u_grid)
]

# Plot
p = plot(u_grid, B_mean,
    ribbon = (B_mean .- B_lower, B_upper .- B_mean),
    fillalpha = 0.2,
    linewidth = 2,
    label = "Average Bernstein CDF (95% bands)",
    xlabel = "u",
    ylabel = "B(u)",
    title = "Average Bernstein CDF with 95% confidence bands")

# Uniform CDF reference
plot!(p, u_grid, u_grid,
    linestyle = :dash,
    linewidth = 2,
    label = "Uniform CDF: B(u) = u")

display(p)