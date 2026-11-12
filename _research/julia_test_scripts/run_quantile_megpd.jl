M = 50
n = 10000
probs = collect(0.001:0.0025:0.999)
kappa_val = 2.0
sigma_val = 1.0
xi_val = 0.5

include("quantile_megpd_markov.jl")

println("Starting simulation:")
println("M = $M, n = $n, kappa = $kappa_val, sigma = $sigma_val, xi = $xi_val")

sim_qq = simulate_megpd_quantiles(
    M, n, probs, kappa_val, sigma_val, xi_val
)

CSV.write("megpd/res/sim_qq_M50_n10k_xi05.csv", Tables.table(sim_qq.Q))