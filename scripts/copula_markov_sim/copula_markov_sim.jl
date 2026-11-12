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

function simulate_series(scenario::AbstractDict, seed::Integer)
    # Extract common parameters
    n = scenario[:n]
    margin_par = [
        scenario[:sigma],
        scenario[:kappa],
        scenario[:xi],
    ]

    # Select the copula-specific parameter
    family = scenario[:family]

    cop_param = scenario[:cop_param]

    # Call the R simulation function
    x = rcopy(
        Vector{Float64},
        R"sim_cm($family, $cop_param, $margin_par, $n, $seed)"
    )

    # Validate the simulation result
    length(x) == n ||
        error("Simulation failed or returned the wrong length")

    all(isfinite, x) ||
        error("Simulation returned non-finite values")

    return x
end


function run_replication(scenario::AbstractDict, rep::Integer)
    # Derive a reproducible seed for this replication
    seed = 10_000 + rep

    # Store scenario information even if the fit fails
    result = Dict{Symbol, Any}(
        :rep => rep,
        :seed => seed,
        :status => "failed",
    )

    # Include all scenario parameters in the result
    merge!(result, scenario)

    try
        # Step 1: simulate
        x = simulate_series(scenario, seed)

        # Step 2: estimate
        fit = MMEGPD.laplace_optimize(
          x, 10;
          spline = :crspline, knot_method = :quantile,
          x_reltol = 1e-5)

        # Step 3: collect estimates
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

# run_replication(all_scenarios[1], 1)

function run_scenario(scenario::AbstractDict)
    nreps = scenario[:MC_reps]
    cop_fam = scenario[:family]
    xi = scenario[:xi]

    println("Running Copula: $cop_fam, xi: $xi \n")
    results = []

    for rep in 1:nreps
        println("  Replication $rep / $nreps")

        result = run_replication(scenario, rep)
        push!(results, result)
    end

    println("Finished Copula: $cop_fam, xi: $xi\n")

    return results
end

# sc_res = run_scenario(all_scenarios[1])

function run_experiment(scenarios)
    all_results = map(run_scenario, scenarios)
    return reduce(vcat, all_results; init=Dict{Symbol, Any}[])
end

results = run_experiment(all_scenarios)

mkpath(datadir())

safesave(
    datadir("copula_markov_sim.jld2"),
    Dict("results" => results)
)