using DataFrames
using Statistics

df = DataFrame(results)

scenario_cols = [
    :family, :cop_param, :n,
    :sigma, :kappa, :xi,
]

function mc_metrics(estimates, truth)
    x = Float64[
        v for v in estimates
        if !ismissing(v) && v isa Real && isfinite(v)
    ]

    if isempty(x)
        return (
            n_valid = 0,
            bias = missing,
            empirical_sd = missing,
            rmse = missing,
        )
    end

    errors = x .- truth

    return (
        n_valid = length(x),
        bias = mean(errors),
        empirical_sd = length(x) > 1 ? std(x) : missing,
        rmse = sqrt(mean(errors .^ 2)),
    )
end

scenario_summary = combine(groupby(df, scenario_cols)) do g
    successful = g[g.status .== "success", :]

    k = mc_metrics(successful.kappa_hat, g.kappa[1])
    s = mc_metrics(successful.sigma_hat, g.sigma[1])
    x = mc_metrics(successful.xi_hat, g.xi[1])

    (
        n_total = nrow(g),
        n_success = nrow(successful),
        failure_rate = 1 - nrow(successful) / nrow(g),

        kappa_n_valid = k.n_valid,
        kappa_bias = k.bias,
        kappa_sd = k.empirical_sd,
        kappa_rmse = k.rmse,

        sigma_n_valid = s.n_valid,
        sigma_bias = s.bias,
        sigma_sd = s.empirical_sd,
        sigma_rmse = s.rmse,

        xi_n_valid = x.n_valid,
        xi_bias = x.bias,
        xi_sd = x.empirical_sd,
        xi_rmse = x.rmse,
    )
end

sort!(scenario_summary, [:family, :xi, :cop_param])

display(scenario_summary)
