using Distributed
using Serialization

if nprocs()==1
    addprocs(49)
end

println("Number of cores: ", nprocs())

@everywhere begin
    using Random
    using MMEGPD
end

# One replicate: simulate + fit. Reseeding per replicate makes this
# deterministic regardless of how @threads schedules tasks.
@everywhere function run_replicate(rep, seed, n, K, kappa, sigma, xi, delta_true; kwargs...)
    Random.seed!(seed + rep)

    Core.println("\n Replicate $rep; sampling data")

    x = MMEGPD.simulate_megpd_chain(n, kappa, sigma, xi, delta_true; show_progress=false)

    Core.println("\n Replicate $rep; optimizing λ")

    la_res = MMEGPD.laplace_optimize(x, K; kwargs...)
    return la_res.best_fit
end

# =============================================================================
# Script entry point
# =============================================================================

if abspath(PROGRAM_FILE) == @__FILE__

    sim_res_path = joinpath(@__DIR__, "res")

    println("Saving directory in $sim_res_path")
    seed=1234
    nrep = nprocs()

    n = 5000
    K = 10
    kappa, sigma, xi = 2.0, 1.0, 0.2
    delta_true = :delta_strong_upper

    @info "Running $nrep replicates on $(nprocs()) processes"

    fits = pmap(
        i -> run_replicate(
            i, seed, n, K, kappa, sigma, xi, delta_true;
            spline=:crspline,
            knot_method=:quantile,
            map=true,
            lower=1.5,
            upper=4.0,
            x_reltol=1e-6
        ),
        1:nrep
    )

    save_file = joinpath(sim_res_path, "cons_par_egpd_crs_quantile_nrep50_n5k_K10_xi02_deltastrongupper.jls")

    Serialization.serialize(save_file, fits)
    @info "Saved results containing the fits object"

end