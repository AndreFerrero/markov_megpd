using DrWatson
@quickactivate "markov_megpd"
using MMEGPD
using RCall
using StatsPlots
using StatsBase
using Random
using JLD2

R"source(\"scripts/winter2016/load_data.R\")"
data = rcopy(R"data")
x = rcopy(R"j.data")

R = x[:,1] .+ x[:,2]
V = x[:,1] ./ x[:,2]

r_grid =
    range(
        minimum(R),
        maximum(R),
        length=500
    )

scatter(
    R,
    log.(V);
    xlabel = "R",
    ylabel = "log(V)",
    label = "Data",
    alpha = 0.4,
    legend = :topright,
)

histogram(R)

K = 10

la_fit = MMEGPD.laplace_optimize(
    x, K;
    verbose=true,
    x_reltol=1e-5,
    kappa_init=30.0,
    sigma_init=10.0
)

@save "data/winter2016/la_fit_mmegpd_K10.jld2" la_fit data x

fit_opt = la_fit.best_fit
fit_opt.theta

MMEGPD.plot_spline(fit_opt)
plot!(title = "λ = $(la_fit.best_lambda)", legend = :topright)
plot!(ylims = (0, 0.3))

fit_opt.basis

delta_fit = delta_spline(fit_opt.basis, fit_opt.beta)

Random.seed!(1234)
M = 150
x0data = rand(x, M)
xrep = MMEGPD.multiple_megpd_chains(M, 2000, fit_opt.theta[1], fit_opt.theta[2], fit_opt.theta[3], delta_fit; x0 = x0data)

@save datadir("winter2016/xrep_la_fit_mmegpd_K10.jld2") xrep x0data
