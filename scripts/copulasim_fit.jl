# scripts/copulasim_fit.jl
#
# Run on the CLUSTER. Loads the simulated data from
#   data/copulasimdata/<experiment_name>/
# fits every replication in parallel, and saves the results in
#   data/copulasimfits/<experiment_name>/
#
# Usage:  julia --project scripts/copulasim_fit.jl <experiment_name>
# No R is needed here.

using DrWatson
@quickactivate "markov_megpd"

include(srcdir("session_info.jl"))
# print_session_info()

println("Loading packages")
flush(stdout)

using MMEGPD
using Dates

include(srcdir("copulasim_funs.jl"))

println("Packages loaded")
flush(stdout)

# if isempty(ARGS)
#     available = isdir(datadir(SIMDATA_SUBDIR)) ? readdir(datadir(SIMDATA_SUBDIR)) : String[]
#     error("Pass the experiment folder name as argument. Available: $available")
# end
# name = ARGS[1]
name = isempty(ARGS) ? first(readdir(datadir(SIMDATA_SUBDIR))) : ARGS[1]

# ---------------------------------------------------------------------------
# Load data (serial), then fit (parallel)
# ---------------------------------------------------------------------------
println("Loading data from ", simdata_dir(name))
flush(stdout)

jobs = load_simulated_data(name)
println("Loaded $(length(jobs)) series")
flush(stdout)

t0 = time()
println("Fit start: ", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
flush(stdout)

results = run_fits(jobs; n_lambda=10, spline=:crspline,
                   knot_method=:quantile, x_reltol=1e-5)

println("Fit finished in $(round(time() - t0, digits=1)) s")
flush(stdout)

dir = save_fits(name, results)
n_ok = count(r -> r["status"] == "success", results)
println("Saved $(length(results)) results ($n_ok successful) in $dir")
flush(stdout)