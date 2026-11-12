# ============================================================
# EGPD: Influence of sigma and kappa
# ============================================================
using DrWatson
@quickactivate "markov_megpd"

using MMEGPD
using Random
using StatsPlots
using Statistics
using Distributions

Random.seed!(1234)

# ------------------------------------------------------------
# 1. General simulation settings
# ------------------------------------------------------------

N = 100_000

# Fixed parameters
sigma_0 = 2.0
xi_0    = 0.1
kappa_0 = 2.0

# Parameter values to investigate
sigma_values = [0.5, 1.0, 2.0, 4.0]

kappa_values = [0.5, 1.0, 2.0, 4.0]

# ------------------------------------------------------------
# 2. EXPERIMENT 1: Influence of sigma
# ------------------------------------------------------------

# Keep kappa and xi fixed, vary sigma.

samples_sigma = Dict(
    s => sample_egpd(N, kappa_0, s, xi_0)
    for s in sigma_values
)

# Evaluation grid
x_grid_sigma = range(
    0.001,
    maximum(maximum.(values(samples_sigma))),
    length=1000
)

# Theoretical densities
p1 = plot(
    xlabel="x",
    ylabel="Density",
    title="Influence of sigma on EGPD density",
    legend=:topright
)

for s in sigma_values

    pdf_values = [
        egpd_lpdf(x, kappa_0, s, xi_0; lpdf=false)
        for x in x_grid_sigma
    ]

    plot!(
        p1,
        x_grid_sigma,
        pdf_values,
        linewidth=2,
        label="sigma = $s"
    )

end

# Theoretical CDFs
p2 = plot(
    xlabel="x",
    ylabel="CDF",
    title="Influence of sigma on EGPD CDF",
    legend=:bottomright
)

for s in sigma_values

    cdf_values = [
        egpd_cdf(x, kappa_0, s, xi_0)
        for x in x_grid_sigma
    ]

    plot!(
        p2,
        x_grid_sigma,
        cdf_values,
        linewidth=2,
        label="sigma = $s"
    )

end

# ------------------------------------------------------------
# 3. EXPERIMENT 2: Influence of kappa on lower-tail behaviour
# ------------------------------------------------------------

# Keep sigma and xi fixed, vary kappa.

samples_kappa = Dict(
    k => sample_egpd(N, k, sigma_0, xi_0)
    for k in kappa_values
)

# Focus on the lower part of the distribution.
# Use the maximum of the selected lower quantiles to define
# a common plotting range.

lower_quantile_max = maximum(
    quantile(samples_kappa[k], 0.10)
    for k in kappa_values
)

x_grid_lower = range(
    1e-6,
    lower_quantile_max,
    length=1000
)

# Theoretical densities near zero
p3 = plot(
    xlabel="x",
    ylabel="Density",
    title="Influence of kappa on lower-tail density",
    legend=:topright
)

for k in kappa_values

    pdf_values = [
        egpd_lpdf(x, k, sigma_0, xi_0; lpdf=false)
        for x in x_grid_lower
    ]

    plot!(
        p3,
        x_grid_lower,
        pdf_values,
        linewidth=2,
        label="kappa = $k"
    )

end

# Theoretical CDFs near zero
p4 = plot(
    xlabel="x",
    ylabel="CDF",
    title="Influence of kappa on lower-tail CDF",
    legend=:bottomright
)

for k in kappa_values

    cdf_values = [
        egpd_cdf(x, k, sigma_0, xi_0)
        for x in x_grid_lower
    ]

    plot!(
        p4,
        x_grid_lower,
        cdf_values,
        linewidth=2,
        label="kappa = $k"
    )

end

# ------------------------------------------------------------
# 4. Empirical lower-tail quantiles
# ------------------------------------------------------------

probs = [0.001, 0.005, 0.01, 0.025, 0.05, 0.10]

println("\nLOWER-TAIL QUANTILES")
println("====================")

println("\nInfluence of sigma:")
println("sigma | " * join(string.(probs), " | "))

for s in sigma_values

    q = quantile.(Ref(samples_sigma[s]), probs)

    println(
        s, " | ",
        join(round.(q, digits=5), " | ")
    )

end

println("\nInfluence of kappa:")
println("kappa | " * join(string.(probs), " | "))

for k in kappa_values

    q = quantile.(Ref(samples_kappa[k]), probs)

    println(
        k, " | ",
        join(round.(q, digits=5), " | ")
    )

end

# ------------------------------------------------------------
# 5. Display plots
# ------------------------------------------------------------

plot(
    p1, p2, p3, p4,
    layout=(2, 2),
    size=(1200, 900)
)

# ============================================================
# END
# ============================================================