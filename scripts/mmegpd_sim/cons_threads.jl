using DrWatson
@quickactivate "markov_megpd"

using Base.Threads
using Random
using JLD2

using MMEGPD

# One replicate: simulate + fit. Reseeding per replicate makes this
# deterministic regardless of how @threads schedules tasks.
function run_replicate(rep, seed, n, K, kappa, xi, delta_true; kwargs...)
    Random.seed!(seed + rep)

    Core.println("\n Replicate $rep; sampling data")
    theta = [kappa, 1.0, xi]
    
    x = MMEGPD.simulate_markov_megpd(n, theta, delta_true; show_progress=false)

    Core.println("\n Replicate $rep; optimizing lambda")

    la_res = MMEGPD.laplace_optimize(x, K; kwargs...)
    return la_res.best_fit
end

# Threaded loop over replicates.
function run_simulation_study(nrep, n, K, kappa, xi, delta_true; seed = 1234, kwargs...)
    fits = Vector{Any}(undef, nrep)
    Threads.@threads for i in 1:nrep
        Core.println("\n Starting replicate $i")
        fits[i] = run_replicate(i, seed, n, K, kappa, xi, delta_true; kwargs...)
        Core.println("\n Finished replicate $i")

    end
    return fits
end

# =============================================================================
# Script entry point
# =============================================================================
println("Saving directory in $(datadir())")

nrep = 100
n = 2000
K = 5
kappa, xi = 2.0, 0.2
delta_true = :delta_exp

@info "Running $nrep replicates on $(Threads.nthreads()) threads"

fits = run_simulation_study(
    nrep, n, K, kappa, xi, delta_true;
    bernstein = true, verbose = false,
    lower = 1.5, upper = 4.0, x_reltol = 1e-6
)

save_file = joinpath(datadir(), "rand_start_cons_bern_egpd_crs_nrep100_n2k_K5_bm5_xi02_deltaexp.jld2")

@save save_file fits nrep n K kappa xi delta_true

@info "Saved results"