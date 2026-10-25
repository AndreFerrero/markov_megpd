# scripts/copulasim_simulate.jl
#
# Run LOCALLY. Simulates data with R (copula / CopulaMarkov via RCall) and
# stores one .jld2 file per (scenario, replication) in
#   data/copulasimdata/<experiment_name>/
# Then copy that folder to the cluster (see the rsync note at the bottom).

using DrWatson
@quickactivate "markov_megpd"

using RCall
using Dates

include(srcdir("copulasim_funs.jl"))

R"""
library(copula)
library(CopulaMarkov)
"""

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

sim_cm <- function(family, cop_param, margin_par, n, seed, max_try = 20) {
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
# Design of the experiment
# ---------------------------------------------------------------------------
n = 2000
mc_reps = 10
true_kappa = 2.0
true_sigma = 1.0
true_xi = [-0.2, 0.2]

families = ["gumbel", "joe", "normal"]
taus = [0.5, 0.75]

common = @strdict n mc_reps true_kappa true_sigma true_xi

all_scenarios = Dict{String,Any}[]
for family in families, tau in taus
    cop_param = rcopy(Float64, R"cop_param_from_tau($family, $tau)")
    append!(all_scenarios, dict_list(merge(
        common,
        Dict("family" => [family], "tau" => [tau], "cop_param" => [cop_param]),
    )))
end

name = experiment_name(; n, mc_reps, families, taus, true_xi)
println("Experiment folder: ", simdata_dir(name))

# ---------------------------------------------------------------------------
# Simulation (R, serial: RCall is not thread-safe)
# ---------------------------------------------------------------------------
function simulate_series(scenario::AbstractDict, seed::Integer)
    n = scenario["n"]
    margin_par = [scenario["true_sigma"], scenario["true_kappa"], scenario["true_xi"]]
    family = scenario["family"]
    cop_param = scenario["cop_param"]

    x = rcopy(Vector{Float64}, R"sim_cm($family, $cop_param, $margin_par, $n, $seed)")

    length(x) == n || error("Simulation failed or returned the wrong length")
    all(isfinite, x) || error("Simulation returned non-finite values")
    return x
end

function simulate_all(name, scenarios)
    t0 = time()
    n_ok, n_fail = 0, 0

    for scenario in scenarios, rep in 1:scenario["mc_reps"]
        meta = replication_meta(scenario, rep)
        try
            x = simulate_series(scenario, meta["seed"])
            save_simulated_data(name, meta, x)
            n_ok += 1
        catch err
            n_fail += 1
            @warn "Simulation failed" family = scenario["family"] tau = scenario["tau"] rep = rep error = sprint(showerror, err)
        end

        if rep % 10 == 0
            println("[$(Dates.format(Dates.now(), "HH:MM:SS"))] ",
                    "$(scenario["family"]), τ=$(scenario["tau"]), rep $rep ",
                    "(ok=$n_ok, failed=$n_fail, $(round(time() - t0, digits=1)) s)")
            flush(stdout)
        end
    end
    return n_ok, n_fail
end

n_ok, n_fail = simulate_all(name, all_scenarios)
println("Done. Saved $n_ok series, $n_fail failed, in ", simdata_dir(name))

# To copy to the cluster (from your local machine):
#   rsync -av data/copulasimdata/<experiment_name> \
#       ferreroa@<cluster>:/home/ferreroa/work_rainstsimu/markov_megpd/data/copulasimdata/