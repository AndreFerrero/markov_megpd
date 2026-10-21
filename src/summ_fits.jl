using Printf

# Sampling distribution of theta = (kappa, sigma, xi) vs. the truth.
function summarize_theta(fits, truth)

    theta = reduce(vcat, (collect(f.theta)' for f in fits))

    mean_est = vec(mean(theta, dims=1))
    sd_est = vec(std(theta, dims=1))
    bias = mean_est .- truth

    ci_lo = [quantile(theta[:, j], 0.025) for j in 1:3]
    ci_hi = [quantile(theta[:, j], 0.975) for j in 1:3]

    return (; mean=mean_est, sd=sd_est, bias, ci_lo, ci_hi)
end


# Fitted delta(r), evaluated on one common r-grid across replicates
# (reuses spline_confidence_band, so it works regardless of whether each
# replicate's basis/knots differ), then averaged pointwise with an MC band.
function summarize_spline(
    fits;
    npoints = 500,
    extrapolate = false,
    level = 0.95
)

    alpha = 1 - level

    # Common grid spanning the union of the spline domains
    r_min = maximum(minimum(f.R) for f in fits)
    r_max = maximum(maximum(f.R) for f in fits)

    r_grid = collect(
        range(r_min, r_max; length = npoints)
    )

    # Evaluate every spline on the SAME grid.
    # Values outside the individual spline domain are missing.
    delta = reduce(
        hcat,
        (
            MMEGPD.spline_confidence_band(
                f,
                r_grid;
                level = level,
                extrapolate = extrapolate
            ).estimate
            for f in fits
        )
    )

    # Pointwise Monte Carlo mean
    delta_mean = [
        mean(skipmissing(view(delta, i, :)))
        for i in axes(delta, 1)
    ]

    # Pointwise Monte Carlo quantiles
    delta_lo = [
        quantile(
            collect(skipmissing(view(delta, i, :))),
            alpha / 2
        )
        for i in axes(delta, 1)
    ]

    delta_hi = [
        quantile(
            collect(skipmissing(view(delta, i, :))),
            1 - alpha / 2
        )
        for i in axes(delta, 1)
    ]

    # Number of simulations contributing at each r
    n_available = [
        count(!ismissing, view(delta, i, :))
        for i in axes(delta, 1)
    ]

    return (
        r_grid = r_grid,
        mean = delta_mean,
        lo = delta_lo,
        hi = delta_hi,
        n_available = n_available
    )
end