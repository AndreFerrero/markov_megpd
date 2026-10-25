# src/copulasim_funs.jl
#
# Reusable, pure-Julia helpers shared by the simulation and fitting scripts.
# No RCall here, so the cluster fitting job never has to start R.
#
# Usage: after `using DrWatson; @quickactivate "markov_megpd"`:
#     include(srcdir("copulasim_funs.jl"))
# The fitting script must also `using MMEGPD` before calling `run_fits`.

using DrWatson
using Dates
using JLD2

const SIMDATA_SUBDIR = "copulasimdata"
const FITS_SUBDIR = "copulasimfits"

# ---------------------------------------------------------------------------
# Naming
# ---------------------------------------------------------------------------

"""Folder name identifying one simulation experiment."""
experiment_name(; n, mc_reps, families, taus, true_xi) =
    "mcreps=$(mc_reps)_n=$(n)_families=$(join(families, "_"))" *
    "_taus=$(join(taus, "_"))_true_xi=$(join(true_xi, "_"))"

"""File name for one (scenario, replication). `meta` must not contain the data."""
simdata_filename(meta::AbstractDict) =
    savename(meta, "jld2"; ignores=("mc_reps", "status", "error"))

simdata_dir(name::AbstractString) = datadir(SIMDATA_SUBDIR, name)
fits_dir(name::AbstractString) = datadir(FITS_SUBDIR, name)

# ---------------------------------------------------------------------------
# Scenarios and metadata
# ---------------------------------------------------------------------------

"""Metadata record for one replication (everything except the data itself)."""
function replication_meta(scenario::AbstractDict, rep::Integer;
                          seed_offset::Integer=10_000)
    meta = Dict{String,Any}("rep" => rep, "seed" => seed_offset + rep)
    merge!(meta, scenario)
    return meta
end

# ---------------------------------------------------------------------------
# Saving / loading simulated data
# ---------------------------------------------------------------------------

"""Save one simulated series with its metadata. Returns the file path."""
function save_simulated_data(name::AbstractString, meta::AbstractDict, x::AbstractVector)
    path = joinpath(simdata_dir(name), simdata_filename(meta))
    d = Dict{String,Any}(meta)
    d["x"] = collect(Float64, x)
    wsave(path, d)
    return path
end

"""Load every simulated series of an experiment. Returns a vector of Dicts
(metadata + "x" + "datafile"), sorted by file name."""
function load_simulated_data(name::AbstractString)
    dir = simdata_dir(name)
    isdir(dir) || error("Data folder not found: $dir")
    files = sort(filter(f -> endswith(f, ".jld2"), readdir(dir; join=true)))
    isempty(files) && error("No .jld2 files in $dir")

    jobs = Dict{String,Any}[]
    for f in files
        d = load(f)                       # Dict{String,Any}
        d["datafile"] = basename(f)
        push!(jobs, d)
    end
    return jobs
end

# ---------------------------------------------------------------------------
# Fitting (thread-safe, pure Julia; requires `using MMEGPD` in the caller)
# ---------------------------------------------------------------------------

"""Fit one replication. `job` contains metadata and "x". Returns a result Dict
(without the data) with estimates and status."""
function fit_replication(job::AbstractDict; n_lambda::Integer=10,
                         spline=:crspline, knot_method=:quantile, x_reltol=1e-5)
    result = Dict{String,Any}(k => v for (k, v) in job if k != "x")
    result["status"] = "failed"

    try
        fit = MMEGPD.laplace_optimize(
            job["x"], n_lambda;
            spline=spline, knot_method=knot_method, x_reltol=x_reltol)

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

"""Fit all jobs in parallel (Threads.@threads) with progress logging.
Returns results in the same order as `jobs`."""
function run_fits(jobs::AbstractVector; kwargs...)
    total_jobs = length(jobs)
    results = Vector{Dict{String,Any}}(undef, total_jobs)

    start = time()
    previous = Ref(start)
    completed = Ref(0)
    lk = ReentrantLock()

    println("Fitting $total_jobs replications on $(Threads.nthreads()) threads")
    flush(stdout)

    Threads.@threads for i in eachindex(jobs)
        result = fit_replication(jobs[i]; kwargs...)
        results[i] = result

        lock(lk) do
            completed[] += 1
            now_t = time()
            ts = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
            println(
                "[$ts] Done: $(result["family"]), τ=$(result["tau"]), " *
                "rep $(result["rep"]) -> $(result["status"]) | " *
                "Completed: $(completed[])/$total_jobs | " *
                "Since previous: $(round(now_t - previous[], digits=2)) s | " *
                "Elapsed: $(round(now_t - start, digits=1)) s"
            )
            flush(stdout)
            previous[] = now_t
        end
    end

    return results
end

"""Save fit results, one file per replication, named like the data files."""
function save_fits(name::AbstractString, results::AbstractVector)
    dir = fits_dir(name)
    for r in results
        wsave(joinpath(dir, r["datafile"]), r)
    end
    return dir
end