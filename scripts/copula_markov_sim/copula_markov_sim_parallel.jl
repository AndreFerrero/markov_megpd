using DrWatson
@quickactivate "markov_megpd"

include(srcdir("session_info.jl"))

print_session_info()

println("Loading packages")
flush(stdout)

using RCall
using MMEGPD

using DrWatson
using RCall

R"""
library(copula)
library(CopulaMarkov)
"""

println("Loaded packages in both languages")
flush(stdout)

R"""
cop_param_from_tau <- function(family, tau) {
    template <- switch(
        family,
        gumbel = gumbelCopula(dim = 2),
        joe    = joeCopula(dim = 2),
        normal = normalCopula(dim = 2)
    )

    round(as.numeric(iTau(template, tau)), 2)
}
"""

n = 2000
mc_reps = 50
true_kappa = 2.0
true_sigma = 1.0
true_xi = [-0.2, 0.2]

common = @strdict n mc_reps true_kappa true_sigma true_xi

families = ["gumbel", "joe", "normal"]
taus = [0.5, 0.75]

all_scenarios = Dict{String,Any}[]

for family in families, tau in taus
    cop_param = rcopy(
        Float64,
        R"cop_param_from_tau($family, $tau)"
    )

    family_scenarios = dict_list(merge(
        common,
        Dict(
            "family" => [family],
            "tau" => [tau],
            "cop_param" => [cop_param],
        )
    ))

    append!(all_scenarios, family_scenarios)
end

for (i, scenario) in enumerate(all_scenarios)
    println("Scenario $i:")
    flush(stdout)
    for (key, value) in scenario
        println("  $key => $value")
        flush(stdout)
    end
    println()
    flush(stdout)
end

R"""
sim_cm <- function(family, cop_param, margin_par, n, seed,
                   max_try = 20) {
  set.seed(seed)
  cop <- switch(family,
    gumbel = copula::gumbelCopula(cop_param),
    joe    = copula::joeCopula(cop_param),
    normal = copula::normalCopula(cop_param, dim = 2))
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

# ---------------------------------------------------------------------------
# Step 1: simulation (R, main thread ONLY -- RCall is not thread-safe)
# ---------------------------------------------------------------------------
function simulate_series(scenario::AbstractDict, seed::Integer)
    n = scenario["n"]
    margin_par = [
        scenario["true_sigma"],
        scenario["true_kappa"],
        scenario["true_xi"],
    ]
    family = scenario["family"]
    cop_param = scenario["cop_param"]

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

# ---------------------------------------------------------------------------
# Step 2 + 3: estimation (pure Julia, thread-safe, run in parallel)
# ---------------------------------------------------------------------------
function fit_replication!(result::Dict{String,Any}, x)
    x === nothing && return result   # simulation failed, keep the failure record

    try
        fit = MMEGPD.laplace_optimize(
            x, 10;
            spline=:crspline, knot_method=:quantile,
            x_reltol=1e-5)

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

# ---------------------------------------------------------------------------
# Experiment: serial simulation, threaded fitting over all (scenario, rep)
# ---------------------------------------------------------------------------
function run_experiment(scenarios)
    experiment_start = time()
    println("Experiment start: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    flush(stdout)

    # Phase 1 (serial): simulate every replication for every scenario
    jobs = Tuple{Dict{String,Any},Any}[]

    for scenario in scenarios
        for rep in 1:scenario["mc_reps"]
            if mod(rep, 10) == 0
                println("Replicate: $rep, Copula: $(scenario["family"]), ξ: $(scenario["true_xi"]), τ: $(scenario["tau"])")
                flush(stdout)
            end
            push!(jobs, simulate_replication(scenario, rep))
        end
    end

    # Phase 2 (parallel): fit. Each job writes only to its own result Dict,
    # so there are no data races. :dynamic balances uneven fit times.
    println("\nFitting $(length(jobs)) replications on $(Threads.nthreads()) threads\n")
    flush(stdout)

    fit_start = time()
    previous_completion = fit_start

    completion_lock = ReentrantLock()

    total_jobs = length(jobs)
    completed_jobs = 0

    Threads.@threads for i in eachindex(jobs)

        result, x = jobs[i]
        fit_replication!(result, x)

        lock(completion_lock) do
            completed_jobs += 1
            remaining_jobs = total_jobs - completed_jobs

            completion_time = time()
            timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
            elapsed_total = completion_time - experiment_start
            elapsed_fit = completion_time - fit_start
            elapsed_since_previous = completion_time - previous_completion

            println(
                "[$timestamp] Done: $(result["family"]), " *
                    "rep $(result["rep"]) -> $(result["status"]) | " *
                    "Completed: $completed_jobs/$total_jobs | " *
                    "Remaining: $remaining_jobs | " *
                    "Since previous fit: $(round(elapsed_since_previous, digits=3)) s | " *
                    "Fitting elapsed: $(round(elapsed_fit, digits=3)) s | " *
                    "Total elapsed: $(round(elapsed_total, digits=3)) s"
            )
            flush(stdout)

            previous_completion = completion_time
        end

    end

    println("Experiment finished: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    flush(stdout)

    # Same ordering as the serial version
    return Dict{String,Any}[first(j) for j in jobs]
end

results = run_experiment(all_scenarios)

ignore_keys = ("mc_reps", "kappa_hat", "sigma_hat", "xi_hat",
    "chosen_lambda", "status")
folder_dir = "copulasim_mcreps=$(mc_reps)_n=$(n)_families=$(join(families, "_"))_taus=$(join(taus, "_"))_true_xi=$(join(true_xi, "_"))"

println("Saving at $folder_dir")

for d in results
    wsave(
        datadir(
            folder_dir,
            savename(d, "jld2"; ignores=ignore_keys)
        ),
        d
    )
end

println("Finished and saved")