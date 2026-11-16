
using DrWatson
@quickactivate "markov_megpd"

using Dates
using RCall
using MMEGPD

include(srcdir("session_info.jl"))

print_session_info()

println("Loading packages")
flush(stdout)

R"""
library(copula)
library(CopulaMarkov)
"""

println("Loaded packages in both languages")
flush(stdout)


# ============================================================
# 1. SCENARIO PARAMETERS
# ============================================================

args = ARGS

if length(args) != 4
    error("""
    Expected 4 arguments:
    family cop_param mc_reps n

    Example:
    julia --threads=4 --project=. scripts/copula_markov_sim_threaded.jl gumbel 0.5 50 2000
    """)
end

# Parameters supplied through the console
family = lowercase(args[1])
cop_param = parse(Float64, args[2])
mc_reps = parse(Int, args[3])
n = parse(Int, args[4])

# Fixed marginal parameters
true_kappa = 2.0
true_sigma = 1.0
true_xi = 0.2

# Bernstein fitting parameters
bernstein = true
bernstein_m = 10

# Other estimation settings
K = 5
spline = :crspline
knot_method = :quantile
x_reltol = 1e-5

# Reproducibility
base_seed = 10_000

# Validate parameters
if !(family in ["gumbel", "joe", "normal"])
    error("Unknown copula family: $family")
end

if mc_reps <= 0 || n <= 0
    error("mc_reps and n must be positive")
end

println("\n========================================")
println("Starting simulation scenario")
println("========================================")
println("Copula family:    $family")
println("Copula parameter: $cop_param")
println("MC repetitions:   $mc_reps")
println("Sample size:      $n")
println("True kappa:       $true_kappa")
println("True sigma:       $true_sigma")
println("True xi:          $true_xi")
println("Bernstein:        $bernstein")
println("Bernstein degree: $bernstein_m")
println("Julia threads:    $(Threads.nthreads())")
println("========================================\n")


# ============================================================
# 3. SCENARIO DICTIONARY
# ============================================================

scenario = Dict{String,Any}(
    "family"       => family,
    "cop_param"    => cop_param,
    "n"            => n,
    "mc_reps"      => mc_reps,
    "true_kappa"   => true_kappa,
    "true_sigma"   => true_sigma,
    "true_xi"      => true_xi,
    "bernstein"    => bernstein,
    "bernstein_m"  => bernstein_m,
    "spline"       => spline,
    "knot_method"  => knot_method,
    "x_reltol"     => x_reltol,
)


# ============================================================
# 4. R SIMULATION FUNCTION
# ============================================================

R"""
sim_cm <- function(family, cop_param, margin_par, n, seed,
                   max_try = 20) {

    set.seed(seed)

    cop <- switch(
        family,
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
# 5. DATA SIMULATION (SERIAL, MAIN THREAD ONLY)
# ============================================================

function simulate_series(scenario::AbstractDict, seed::Integer)

    n = scenario["n"]

    margin_par = [
        scenario["true_sigma"],
        scenario["true_kappa"],
        scenario["true_xi"],
    ]

    family = scenario["family"]
    cop_param = scenario["cop_param"]

    # RCall: never execute this function from a worker thread
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


# Simulate one replication.
# Returns (result, x), with x === nothing on failure.

function simulate_replication(scenario::AbstractDict, rep::Integer)

    seed = 10_000 + rep

    result = Dict{String,Any}(
        "rep" => rep,
        "seed" => seed,
        "status" => "failed",
    )

    merge!(result, scenario)

    x = nothing

    try
        x = simulate_series(scenario, seed)

    catch err
        result["error"] = sprint(showerror, err)
    end

    return result, x
end


# ============================================================
# 6. MODEL FITTING (PURE JULIA, MULTITHREADED)
# ============================================================

function fit_replication!(
    result::Dict{String,Any},
    x
)

    # Simulation failed: retain failure record
    x === nothing && return result

    try

        fit = MMEGPD.laplace_optimize(
            x, K;
            bernstein = result["bernstein"],
            bernstein_m = result["bernstein_m"],
            spline = result["spline"],
            knot_method = result["knot_method"],
            x_reltol = result["x_reltol"]
        )

        best_fit = fit.best_fit

        result["fit_obj"] = best_fit

        theta_hat = best_fit.theta

        result["kappa_hat"] = theta_hat[1]
        result["sigma_hat"] = theta_hat[2]
        result["xi_hat"] = theta_hat[3]

        result["chosen_lambda"] = fit.best_lambda
        result["status"] = "success"

    catch err

        result["error"] = sprint(showerror, err)

    end

    return result
end


# ============================================================
# 7. RUN THE SCENARIO
# ============================================================

function run_scenario(scenario::AbstractDict)

    experiment_start = time()

    println("Experiment start: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    flush(stdout)

    # --------------------------------------------------------
    # Phase 1: SERIAL DATA SIMULATION
    # --------------------------------------------------------

    println("\nPhase 1: Simulating datasets (serial, R)")
    flush(stdout)

    jobs = Tuple{Dict{String,Any},Any}[]

    for rep in 1:scenario["mc_reps"]

        if mod(rep, 10) == 0 || rep == 1
            println(
                "Simulating replication $rep / $(scenario["mc_reps"])"
            )
            flush(stdout)
        end

        push!(
            jobs,
            simulate_replication(scenario, rep)
        )
    end

    println("\nSimulation completed: $(length(jobs)) datasets")
    flush(stdout)


    # --------------------------------------------------------
    # Phase 2: PARALLEL MODEL FITTING (PURE JULIA)
    # --------------------------------------------------------

    println(
        "\nPhase 2: Fitting $(length(jobs)) replications " *
        "using $(Threads.nthreads()) Julia threads\n"
    )
    flush(stdout)

    fit_start = time()
    previous_completion = fit_start

    completion_lock = ReentrantLock()

    total_jobs = length(jobs)
    completed_jobs = 0

    Threads.@threads for i in eachindex(jobs)

        result, x = jobs[i]

        fit_replication!(result, x)

        # Synchronise progress reporting only.
        # Each thread writes to its own result dictionary.

        lock(completion_lock) do

            completed_jobs += 1

            remaining_jobs = total_jobs - completed_jobs

            completion_time = time()

            timestamp = Dates.format(
                now(),
                "yyyy-mm-dd HH:MM:SS.sss"
            )

            elapsed_total = completion_time - experiment_start
            elapsed_fit = completion_time - fit_start

            elapsed_since_previous =
                completion_time - previous_completion

            println(
                "[$timestamp] Done: $(result["family"]), " *
                "rep $(result["rep"]) -> $(result["status"]) | " *
                "Completed: $completed_jobs/$total_jobs | " *
                "Remaining: $remaining_jobs | " *
                "Since previous fit: " *
                "$(round(elapsed_since_previous, digits=3)) s | " *
                "Fitting elapsed: $(round(elapsed_fit, digits=3)) s | " *
                "Total elapsed: $(round(elapsed_total, digits=3)) s"
            )

            flush(stdout)

            previous_completion = completion_time
        end
    end

    # --------------------------------------------------------
    # Phase 3: COMPLETION SUMMARY
    # --------------------------------------------------------

    successful = count(
        result -> result["status"] == "success",
        first.(jobs)
    )

    println("\n========================================")
    println("Scenario completed")
    println("========================================")
    println("Copula: $family")
    println("Successful fits: $successful / $total_jobs")
    println("Failed fits: $(total_jobs - successful)")
    println(
        "Total elapsed: " *
        "$(round(time() - experiment_start, digits=2)) seconds"
    )
    println("========================================\n")

    flush(stdout)

    # Preserve original replication ordering
    return Dict{String,Any}[first(j) for j in jobs]
end


# ============================================================
# 8. EXECUTE SCENARIO
# ============================================================

results = run_scenario(scenario)


# ============================================================
# 9. SAVE RESULTS
# ============================================================

folder_dir = "copulasim" *
    "_family=$(family)" *
    "_cop_param=$(cop_param)" *
    "_MC=$(mc_reps)" *
    "_n=$(n)" *
    "_xi=$(join(true_xi, "_"))" *
    "_bernstein=$(bernstein_m)"

output_dir = datadir(folder_dir)

mkpath(output_dir)

println("Saving results at: $output_dir")
flush(stdout)

# Save each replication separately
ignore_keys = (
    "fit_obj",
    "error",
    "kappa_hat",
    "sigma_hat",
    "xi_hat",
    "chosen_lambda",
    "status",
)

for d in results

    wsave(
        datadir(
            folder_dir,
            savename(d, "jld2"; ignores=ignore_keys)
        ),
        d
    )
end

println("\nFinished and saved all replications.")
flush(stdout)
