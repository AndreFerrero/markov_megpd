
using Base.Threads
using Random
using Serialization

using MMEGPD

# One replicate: simulate + fit. Reseeding per replicate makes this
# deterministic regardless of how @threads schedules tasks.
function run_replicate(rep, seed, n, K, kappa, sigma, xi, delta_true; kwargs...)
    Random.seed!(seed + rep)

    Core.println("\n Replicate $rep; sampling data")

    x = MMEGPD.simulate_megpd_chain(n, kappa, sigma, xi, delta_true; show_progress=false)

    Core.println("\n Replicate $rep; optimizing λ")

    la_res = MMEGPD.laplace_optimize(x, K; kwargs...)
    return la_res.best_fit
end

# Threaded loop over replicates.
function run_simulation_study(nrep, n, K, kappa, sigma, xi, delta_true; seed = 1234, kwargs...)
    fits = Vector{Any}(undef, nrep)
    Threads.@threads for i in 1:nrep
        Core.println("\n Starting replicate $i")
        fits[i] = run_replicate(i, seed, n, K, kappa, sigma, xi, delta_true; kwargs...)
        Core.println("\n Finished replicate $i, remaining: $(nrep-i)/$nrep")

    end
    return fits
end

# =============================================================================
# Script entry point
# =============================================================================

if abspath(PROGRAM_FILE) == @__FILE__

    sim_res_path = joinpath(@__DIR__, "res")

    println("Saving directory in $sim_res_path")
    nrep = 50
    n = 5000
    K = 10
    kappa, sigma, xi = 2.0, 1.0, 0.5
    delta_true = :delta_strong_upper

    @info "Running $nrep replicates on $(Threads.nthreads()) threads"

    fits = run_simulation_study(
        nrep, n, K, kappa, sigma, xi, delta_true;
        spline = :crspline, knot_method = :quantile,
        map = true, lower = 1.5, upper = 4.0, x_reltol = 1e-6
    )

    save_file = joinpath(sim_res_path, "cons_egpd_crs_nrep50_n5k_K10_xi05_deltastrongupper.jls")

    Serialization.serialize(save_file, fits)
    @info "Saved results"

end