# ============================================================================
# run_scenario.jl : defines ONE scenario and runs it with sim_core.jl
#
#   julia run_scenario.jl                      # pilot, default scenario
#   julia run_scenario.jl joe 4.0 0.2          # family param xi  (pilot)
#   SIM_MODE=full julia run_scenario.jl joe 4.0 0.2
#   SIM_MODE=full julia --threads=4 run_scenario.jl normal 0.8 -0.3
#
# The family/param/xi can be passed as arguments, which makes it easy to
# launch the whole grid from a shell loop or a cluster job array.
# ============================================================================
include("sim_core.jl")

mode   = get(ENV, "SIM_MODE", "pilot")
family = length(ARGS) >= 1 ? ARGS[1] : "gumbel"
param  = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 2.0
xi     = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.2

# all fitting variants; ASSUMPTION: keyword names beyond the first must be adapted to MMEGPD
all_configs = Any[
    (label = "cr_quantile", kw = (spline = :crspline, knot_method = :quantile)),
    (label = "cr_even",     kw = (spline = :crspline, knot_method = :even)),
    (label = "ps1",         kw = (spline = :pspline,  penalty_order = 1)),
    (label = "ps2",         kw = (spline = :pspline,  penalty_order = 2)),
]

sc = if mode == "pilot"
    Scenario(
        family = family, param = param, xi = xi,
        n = 1000, nrep = 5,
        configs = all_configs[1:1],
        outdir = "sim_pilot",
        keep_fits = 5,
    )
elseif mode == "full"
    Scenario(
        family = family, param = param, xi = xi,
        n = 2000, nrep = 200,
        configs = all_configs,            # compare all spline variants
        outdir = "sim_results",
        keep_fits = 10,
        use_threads = Threads.nthreads() > 1,
    )
else
    error("SIM_MODE must be 'pilot' or 'full', got '$mode'")
end

df = run_scenario(sc; overwrite = (mode == "pilot"))
println(summarize_results(sc.outdir))