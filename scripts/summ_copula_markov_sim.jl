using DataFrames
using Statistics
using DrWatson

include(srcdir("copulasim_funs.jl"))

name = isempty(ARGS) ? first(readdir(datadir(FITS_SUBDIR))) : ARGS[1]

# df = collect_results(datadir("copulasim"))
load_simulated_data(name)

df = collect_results(datadir("copulasimfits/$name"))

df[!,Not([:fit_obj, :status, :seed])]

xi_results = select(
    df,
    :family,
    :tau,
    :true_xi,
    :xi_hat,
    :rep
)

xi_summary = combine(
    groupby(xi_results, [:family, :tau, :true_xi])
) do g
    estimates = g.xi_hat
    truth = first(g.true_xi)

    DataFrame(
        n_reps = length(estimates),
        mean_xi_hat = mean(estimates),
        bias = mean(estimates) - truth,
        sd_xi_hat = std(estimates),
        rmse = sqrt(mean((estimates .- truth).^2))
    )
end

sort!(xi_summary, [:true_xi, :tau, :family])

display(xi_summary)
