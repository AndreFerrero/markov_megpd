using RCall
using MMEGPD

all_scenarios = dict_list(Dict(
    :n => [2000], # Sample size
    :MC_reps => [1000], # Number of MC simulations per scenario,
    :kappa => [2.0],
    :sigma => [1.0],
    :xi => [-0.3, 0.2, 0.5],
    
))


R"""
library(CopulaMarkov)
"""

R"""
sim_cm <- function(family, param, margin_par, n, seed, max_try = 20) {
  set.seed(seed)
  cop <- switch(family,
    gumbel = copula::gumbelCopula(param),
    joe    = copula::joeCopula(param),
    normal = copula::normalCopula(param, dim = 2))
  mod <- copula_markov(egpd_margin(margin_par), cop)
  for (i in seq_len(max_try)) {
    out <- try(simulate(mod, nsim = n), silent = TRUE)
    if (!inherits(out, "try-error") && all(is.finite(out$X))) return(as.numeric(out$X))
  }
  numeric(0)
}
"""

function simulate_series(scenario, seed::Int)
    mp = [sc.sigma, sc.kappa, sc.xi]
    fam, par, n = sc.family, sc.param, sc.n
    return rcopy(Vector{Float64}, R"sim_cm($fam, $par, $mp, $n, $seed)")
end
