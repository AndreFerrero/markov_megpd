using DrWatson
@quickactivate "markov_megpd"

using RCall
using MMEGPD

# ============================================================
# 1. SCENARIO PARAMETERS
# ============================================================

args = ARGS

if length(args) != 4
    error("""
    Expected 4 arguments:
    family cop_param MC_reps n

    Example:
    julia scripts/copula_markov_sim.jl gumbel 2.0 100 100000
    """)
end

# Parameters supplied through the console
family    = lowercase(args[1])
cop_param = parse(Float64, args[2])
MC_reps   = parse(Int, args[3])
n         = parse(Int, args[4])

# Fixed marginal parameters
kappa = 2.0
sigma = 1.0
xi    = 0.2

# Bernstein fitting parameters
bernstein   = true
bernstein_m = 10

# Other estimation settings
K = 10
spline      = :crspline
knot_method = :quantile
x_reltol    = 1e-5

# Validate parameters
if !(family in ["gumbel", "joe", "normal"])
    error("Unknown copula family: $family")
end

if MC_reps <= 0 || n <= 0
    error("MC_reps and n must be positive")
end

if family == "gumbel" && cop_param < 1
    error("Gumbel copula requires cop_param >= 1")
elseif family == "joe" && cop_param < 1
    error("Joe copula requires cop_param >= 1")
elseif family == "normal" && abs(cop_param) >= 1
    error("Normal copula requires -1 < cop_param < 1")
end

# Scenario identification
scenario_name = "$(family)$(cop_param)_MC$(MC_reps)_n$(n)_bernstein$(bernstein_m)"

println("\n========================================")
println("Starting simulation scenario")
println("========================================")
println("Copula family:    $family")
println("Copula parameter: $cop_param")
println("MC repetitions:   $MC_reps")
println("Sample size:      $n")
println("kappa:            $kappa")
println("sigma:            $sigma")
println("xi:               $xi")
println("Bernstein:        $bernstein")
println("Bernstein degree: $bernstein_m")
println("Scenario:         $scenario_name")
println("========================================\n")


# ============================================================
# 2. SCENARIO DICTIONARY
# ============================================================

scenario = Dict{Symbol, Any}(
    :family       => family,
    :cop_param    => cop_param,
    :MC_reps      => MC_reps,
    :n            => n,
    :kappa        => kappa,
    :sigma        => sigma,
    :xi           => xi,
    :bernstein    => bernstein,
    :bernstein_m  => bernstein_m,
    :spline       => spline,
    :knot_method  => knot_method,
    :x_reltol     => x_reltol,
)


# ============================================================
# 3. R SIMULATION FUNCTION
# ============================================================

R"""
library(CopulaMarkov)

sim_cm <- function(family, cop_param, margin_par, n, seed, max_try = 20) {

  set.seed(seed)

  cop <- switch(family,
    gumbel = copula::gumbelCopula(cop_param),
    joe    = copula::joeCopula(cop_param),
    normal = copula::normalCopula(cop_param, dim = 2)
  )

  mod <- copula_markov(egpd_margin(margin_par), cop)

  for (i in seq_len(max_try)) {

    out <- try(simulate(mod, nsim = n), silent = TRUE)

    if (!inherits(out, "try-error") && all(is.finite(out$X))) {
      return(as.numeric(out$X))
    }
  }

  numeric(0)
}
"""


# ============================================================
# 4. DATA SIMULATION
# ============================================================

function simulate_series(scenario::AbstractDict, seed::Integer)

    n = scenario[:n]

    margin_par = [
        scenario[:sigma],
        scenario[:kappa],
        scenario[:xi],
    ]

    family = scenario[:family]
    cop_param = scenario[:cop_param]

    # Call the R simulation function
    x = rcopy(
        Vector{Float64},
        R"sim_cm($family, $cop_param, $margin_par, $n, $seed)"
    )

    # Validate simulation
    length(x) == n ||
        error("Simulation failed or returned the wrong length")

    all(isfinite, x) ||
        error("Simulation returned non-finite values")

    return x
end


# ============================================================
# 5. ONE MONTE CARLO REPLICATION
# ============================================================

function run_replication(scenario::AbstractDict, rep::Integer)

    # Reproducible seed
    seed = 10_000 + rep

    # Store scenario information even if fitting fails
    result = Dict{Symbol, Any}(
        :rep => rep,
        :seed => seed,
        :status => "failed",
    )

    merge!(result, scenario)

    try

        # Step 1: simulate data
        x = simulate_series(scenario, seed)

        # Step 2: estimate MEGPD model
        fit = MMEGPD.laplace_optimize(
            x, K;
            bernstein = scenario[:bernstein],
            bernstein_m = scenario[:bernstein_m],
            spline = scenario[:spline],
            knot_method = scenario[:knot_method],
            x_reltol = scenario[:x_reltol]
        )

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


# ============================================================
# 6. RUN THE SCENARIO
# ============================================================

function run_scenario(scenario::AbstractDict)

    nreps = scenario[:MC_reps]
    cop_fam = scenario[:family]
    cop_param = scenario[:cop_param]
    xi = scenario[:xi]

    println("Running Copula: $cop_fam")
    println("Copula parameter: $cop_param")
    println("xi: $xi\n")

    results = Dict{Symbol, Any}[]

    for rep in 1:nreps

        println("  Replication $rep / $nreps")

        result = run_replication(scenario, rep)

        push!(results, result)

        if result[:status] == "failed"
            println("  Failed: ", result[:error])
        end
    end

    successful = count(r -> r[:status] == "success", results)

    println("\nFinished scenario: $cop_fam")
    println("Successful fits: $successful / $nreps\n")

    return results
end


# ============================================================
# 7. EXECUTE AND SAVE
# ============================================================

results = run_scenario(scenario)

mkpath(datadir())

output_file = datadir("copula_markov_$(scenario_name).jld2")

safesave(
    output_file,
    Dict(
        "results" => results,
        "scenario" => scenario,
    )
)

println("Results saved to: $output_file")
