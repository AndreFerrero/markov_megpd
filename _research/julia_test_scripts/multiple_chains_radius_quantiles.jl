
using MMEGPD

using Random
using Statistics
using StatsBase
using Plots
using Printf

Random.seed!(1234)

# ============================================================
# PARAMETERS
# ============================================================

# Number of independent chains
n_chains = 50

# Number of iterations per chain
n_steps = 5000

# MMEGPD parameters
kappa = 2.0
sigma = 1.0
xi = 0.5

# Delta function
delta_name = :delta_strong_upper

# Initial value
x0 = 1.0

# Fraction of each chain discarded as burn-in
burn_in_prop = 0.0


# ============================================================
# QUANTILE PROBABILITIES
# ============================================================

probs = [
    0.05,
    0.10,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.99,
    0.999
]


# ============================================================
# SIMULATE CHAINS
# ============================================================

println("Simulating $n_chains chains...")
println("n_steps = $n_steps")
println("burn-in = $(100 * burn_in_prop)%")
println()

chains = MMEGPD.multiple_megpd_chains(
    n_chains,
    n_steps,
    kappa,
    sigma,
    xi,
    delta_name;
    x0=x0,
    burn_in_prop=burn_in_prop
)

println("Simulation finished.")
println()


# ============================================================
# CONSTRUCT THE RADIUS
# ============================================================

# R_t = X_t + X_{t-1}
#
# Each column of chains is one simulated chain.
#
# The radius therefore has n_steps - 1 observations
# per chain when burn_in_prop = 0.

radius_quantiles = Matrix{Float64}(
    undef,
    length(probs),
    n_chains
)

maxR = Vector(undef, n_chains)

for j in 1:n_chains

    x = chains[:, j]

    # Radius
    R = x[2:end] .+ x[1:end-1]

    println("Chain $j maximum:", maximum(R))
    maxR[j] = maximum(R)

    # Quantiles of the radius
    radius_quantiles[:, j] = quantile(
        R,
        probs
    )
end

# Checking radial maxima distribution (iid vs megpd)

# histogram(maxR)

# maxRltx = sum(maxR .< 200)
# maxRltx / n_chains

# MMEGPD.egpd_max_cdf(200, n_steps, kappa, sigma, xi)

# maxRiid = Vector(undef, n_chains)
# for j in 1:n_chains
#     Riid = MMEGPD.sample_egpd(5000, 2.0, 1.0, 0.5)
#     maxRiid[j] = maximum(Riid)
# end

# maxRiidltx = sum(maxRiid .< 200)
# maxRiidltx / n_chains

# ============================================================
# AVERAGE RADIUS QUANTILES ACROSS CHAINS
# ============================================================

mean_radius_quantiles = vec(
    mean(radius_quantiles, dims=2)
)


# ============================================================
# EXACT EGPD QUANTILES
# ============================================================

theoretical_quantiles = [
    MMEGPD.egpd_quantile(
        p,
        kappa,
        sigma,
        xi
    )
    for p in probs
]


# ============================================================
# COMPARISON TABLE
# ============================================================

println()
println("==============================================================")
println("EGPD vs MEGPD RADIUS QUANTILES")
println("==============================================================")

@printf(
    "%8s %18s %18s %18s %18s\n",
    "p",
    "EGPD",
    "Mean radius",
    "Difference",
    "Relative"
)

println("--------------------------------------------------------------")

for i in eachindex(probs)

    q_theory = theoretical_quantiles[i]
    q_radius = mean_radius_quantiles[i]

    difference = q_radius - q_theory
    relative = difference / q_theory

    @printf(
        "%8.3f %18.6f %18.6f %18.6f %18.4f\n",
        probs[i],
        q_theory,
        q_radius,
        difference,
        relative
    )
end

println("==============================================================")


# ============================================================
# Q-Q PLOT: EGPD vs RADIUS
# ============================================================

p = scatter(
    theoretical_quantiles,
    mean_radius_quantiles;
    xlabel = "Theoretical EGPD quantiles",
    ylabel = "Mean empirical radius quantiles",
    markersize = 5,
    legend = false
)

# 45-degree reference line
q_min = min(
    minimum(theoretical_quantiles),
    minimum(mean_radius_quantiles)
)

q_max = max(
    maximum(theoretical_quantiles),
    maximum(mean_radius_quantiles)
)

plot!(
    p,
    [q_min, q_max],
    [q_min, q_max];
    linestyle = :dash,
    linewidth = 2
)

display(p)
