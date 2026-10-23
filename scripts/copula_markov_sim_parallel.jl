using DrWatson
@quickactivate "markov_megpd"

using RCall
using MMEGPD

common = Dict(
    :n       => [100],
    :MC_reps => [1],
    :kappa   => [2.0],
    :sigma   => [1.0],
    :xi      => [0.2],
)

# Gumbel and Joe copulas
scenarios_archimedean = dict_list(merge(
    common,
    Dict(
        :family    => ["gumbel"],
        :cop_param => [2.0],
    )
))

# Gaussian copula
scenarios_normal = dict_list(merge(
    common,
    Dict(
        :family    => ["normal"],
        :cop_param => [0.5],
    )
))

# Combine all scenarios
all_scenarios = vcat(
    scenarios_archimedean,
    scenarios_normal,
)

R"""
library(CopulaMarkov)
"""

R"""
sim_cm <- function(family, cop_param, margin_par, n, seed, max_try = 20) {
  set.seed(seed)
  cop <- switch(family,
    gumbel = copula::gumbelCopula(cop_param),
    joe    = copula::joeCopula(cop_param),
    normal = copula::normalCopula(cop_param, dim = 2))
  mod <- copula_markov(egpd_margin(margin_par), cop)
  for (i in seq_len(max_try)) {
    out <- try(simulate(mod, nsim = n), silent = TRUE)
    if (!inherits(out, "try-error") && all(is.finite(out$X))) return(as.numeric(out$X))
  }
  numeric(0)
}
"""

# ---------------------------------------------------------------------------
# Step 1: simulation (R, main thread ONLY -- RCall is not thread-safe)
# ---------------------------------------------------------------------------
function simulate_series(scenario::AbstractDict, seed::Integer)
    n = scenario[:n]
    margin_par = [
        scenario[:sigma],
        scenario[:kappa],
        scenario[:xi],
    ]
    family    = scenario[:family]
    cop_param = scenario[:cop_param]

    x = rcopy(
        Vector{Float64},
        R"sim_cm($family, $cop_param, $margin_par, $n, $seed)"
    )

    length(x) == n ||
        error("Simulation failed or returned the wrong length")

    all(isfinite, x) ||
        error("Simulation returned non-finite values")

    return x
end

# Simulate one replication. Returns (result, x), with x === nothing on failure.
function simulate_replication(scenario::AbstractDict, rep::Integer)
    seed = 10_000 + rep

    result = Dict{Symbol, Any}(
        :rep => rep,
        :seed => seed,
        :status => "failed",
    )
    merge!(result, scenario)

    x = nothing
    try
        x = simulate_series(scenario, seed)
    catch err
        result[:error] = sprint(showerror, err)
    end

    return result, x
end

# ---------------------------------------------------------------------------
# Step 2 + 3: estimation (pure Julia, thread-safe, run in parallel)
# ---------------------------------------------------------------------------
function fit_replication!(result::Dict{Symbol, Any}, x)
    x === nothing && return result   # simulation failed, keep the failure record

    try
        fit = MMEGPD.laplace_optimize(
          x, 10;
          spline = :crspline, knot_method = :quantile,
          x_reltol = 1e-5)

        best_fit = fit.best_fit
        result[:fit_obj] = best_fit

        theta_hat = best_fit.theta
        result[:kappa_hat] = theta_hat[1]
        result[:sigma_hat] = theta_hat[2]
        result[:xi_hat]    = theta_hat[3]

        result[:chosen_lambda] = fit.best_lambda
        result[:status] = "success"
    catch err
        result[:error] = sprint(showerror, err)
    end

    return result
end

# ---------------------------------------------------------------------------
# Experiment: serial simulation, threaded fitting over all (scenario, rep)
# ---------------------------------------------------------------------------
function run_experiment(scenarios)
    # Phase 1 (serial): simulate every replication of every scenario
    jobs = Tuple{Dict{Symbol, Any}, Any}[]
    for scenario in scenarios
        println("Simulating Copula: $(scenario[:family]), xi: $(scenario[:xi])")
        
        for rep in 1:scenario[:MC_reps]
            push!(jobs, simulate_replication(scenario, rep))
        end
    end

    # Phase 2 (parallel): fit. Each job writes only to its own result Dict,
    # so there are no data races. :dynamic balances uneven fit times.
    println("\nFitting $(length(jobs)) replications on $(Threads.nthreads()) threads\n")
    Threads.@threads :dynamic for i in eachindex(jobs)
        result, x = jobs[i]
        fit_replication!(result, x)
        println("  Done: $(result[:family]), rep $(result[:rep]) -> $(result[:status])")
    end

    # Same ordering as the serial version
    return Dict{Symbol, Any}[first(j) for j in jobs]
end

results = run_experiment(all_scenarios)

mkpath(datadir())

safesave(
    datadir("copula_markov_sim.jld2"),
    Dict("results" => results)
)