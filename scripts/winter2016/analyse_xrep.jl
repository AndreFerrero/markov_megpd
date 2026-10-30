using DrWatson
@quickactivate "markov_megpd"

using JLD2
using RCall
using Statistics
using Plots

include(srcdir("inter_theta_clust.jl"))   # intervals_estimator, decluster_intervals

@load datadir("winter2016/la_fit_mmegpd_K10.jld2") la_fit data
@load datadir("winter2016/xrep_la_fit_mmegpd_K10.jld2") xrep x0data

M = size(xrep, 2)                          # number of replicates (columns of xrep)

############################################################
## HELPERS
############################################################

# NaN-aware quantile (replicates with no exceedances give NaN statistics)
function qnan(v, p)
    w = filter(!isnan, v)
    return isempty(w) ? NaN : quantile(w, p)
end

# NaN-aware mean (NaN if every replicate is NaN)
function mnan(v)
    w = filter(!isnan, v)
    return isempty(w) ? NaN : mean(w)
end

# Row-wise summary of a (levels x replicates) matrix
function rep_summary(rep::AbstractMatrix; lo = 0.05, hi = 0.95)
    rows = axes(rep, 1)
    return (
        mean   = [mnan(rep[i, :])          for i in rows],
        median = [qnan(rep[i, :], 0.5)     for i in rows],
        lower  = [qnan(rep[i, :], lo)      for i in rows],
        upper  = [qnan(rep[i, :], hi)      for i in rows],
        nvalid = [count(!isnan, rep[i, :]) for i in rows],   # replicates used per level
    )
end

# Observed vs simulated median and 5-95% band, one panel
function band_plot(levels, obs, rep; xlabel = "", ylabel = "", title = "",
                   legend = false)
    s = rep_summary(rep)
    p = plot(levels, obs;
        marker = :circle, linewidth = 2, label = "Observed",
        xlabel = xlabel, ylabel = ylabel, title = title, legend = legend)
    plot!(p, levels, s.mean; marker = :utriangle, linewidth = 2,
          linestyle = :dot, label = "Simulated mean")
    plot!(p, levels, s.lower; linestyle = :dash, label = "Simulated 5th percentile")
    plot!(p, levels, s.upper; linestyle = :dash, label = "Simulated 95th percentile")
    return p
end

# Observed value in the replicate interval?
inband(obs, s, i) = s.lower[i] <= obs[i] <= s.upper[i]

############################################################
## QUANTILES
############################################################

probs = [0.90, 0.95, 0.975, 0.99, 0.995, 0.999]

q_obs = quantile(data, probs)
q_rep = reduce(hcat, [quantile(view(xrep, :, j), probs) for j in 1:M])  # probs x M
q_sum = rep_summary(q_rep)

for i in eachindex(probs)
    println(
        "p = ", probs[i],
        " | observed = ", round(q_obs[i], digits = 4),
        " | simulated mean = ", round(q_sum.mean[i], digits = 4),
        " | simulated 90% interval = [",
        round(q_sum.lower[i], digits = 4), ", ",
        round(q_sum.upper[i], digits = 4), "]"
    )
end

p_quantiles = band_plot(probs, q_obs, q_rep;
    xlabel = "Quantile probability", ylabel = "Quantile value",
    title = "Upper-tail quantile reproduction", legend = :topleft)
display(p_quantiles)

############################################################
## EXCEEDANCE FREQUENCY
############################################################

# Thresholds are fixed from the OBSERVED data and applied to every replicate
threshold_probs = [0.95, 0.99, 0.995]
thresholds = quantile(data, threshold_probs)

exc_obs = [mean(data .> u) for u in thresholds]
exc_rep = [mean(view(xrep, :, j) .> u) for u in thresholds, j in 1:M]  # thresholds x M

p_exceedance = band_plot(threshold_probs, exc_obs, exc_rep;
    xlabel = "Threshold quantile", ylabel = "Fraction exceeding threshold",
    title = "Exceedance frequency", legend = :topright)
display(p_exceedance)

############################################################
## CLUSTERS OF EXTREMES (interexceedance times, Ferro & Segers)
############################################################

# Cluster summary of one series at threshold u
function cluster_stats(x, u)
    if !any(>(u), x)
        return (nexc = 0.0, theta = NaN, nclusters = 0.0, mean_size = NaN,
                max_size = NaN, mean_duration = NaN, run_length = NaN)
    end
    res       = decluster_intervals(x, u)
    sizes     = length.(res.clusters)                       # exceedances per cluster
    durations = [c[end] - c[1] + 1 for c in res.clusters]   # time span of each cluster
    return (
        nexc          = float(sum(sizes)),
        theta         = res.theta,
        nclusters     = float(res.nclusters),
        mean_size     = mean(sizes),
        max_size      = float(maximum(sizes)),
        mean_duration = mean(durations),
        run_length    = res.run_length,
    )
end

obs_cl = [cluster_stats(data, u) for u in thresholds]
rep_cl = [cluster_stats(view(xrep, :, j), u) for u in thresholds, j in 1:M]

# Extract one statistic as vector (observed) / matrix (thresholds x M, replicates)
obs_stat(f) = [getfield(s, f) for s in obs_cl]
rep_stat(f) = [getfield(s, f) for s in rep_cl]

cluster_labels = [
    (:theta,         "Extremal index (intervals estimator)"),
    (:nclusters,     "Number of clusters"),
    (:mean_size,     "Mean cluster size (exceedances)"),
    (:max_size,      "Maximum cluster size"),
    (:mean_duration, "Mean cluster duration (time steps)"),
]

# Console summary
for (f, label) in cluster_labels
    o = obs_stat(f)
    s = rep_summary(rep_stat(f))
    println("\n", label)
    for i in eachindex(threshold_probs)
        println(
            "  threshold q", threshold_probs[i],
            " | observed = ", round(o[i], digits = 3),
            " | simulated median = ", round(s.median[i], digits = 3),
            " | simulated mean = ", round(s.mean[i], digits = 3),
            " | 90% interval = [", round(s.lower[i], digits = 3), ", ",
            round(s.upper[i], digits = 3), "]",
            " | valid reps = ", s.nvalid[i],
            inband(o, s, i) ? "" : "  <-- outside"
        )
    end
end

# One panel per statistic
cluster_panels = [
    band_plot(threshold_probs, obs_stat(f), rep_stat(f);
        xlabel = "Threshold quantile", ylabel = label, title = label,
        legend = (f == :theta ? :best : false))
    for (f, label) in cluster_labels
]
p_clusters = plot(cluster_panels...; layout = (3, 2), size = (1000, 1100))
display(p_clusters)

############################################################
## EXTREMOGRAM
############################################################
# rho_u(h) = P(X_{t+h} > u | X_t > u), estimated empirically for lags h = 1..H
# at the same thresholds as above (fixed from the observed data).

function extremogram(x, u, lags)
    ex = x .> u
    n  = length(x)
    return [
        begin
            den = count(view(ex, 1:n-h))                   # exceedances at time t
            num = count(i -> ex[i] && ex[i+h], 1:n-h)      # ... also exceeding at t+h
            den == 0 ? NaN : num / den
        end
        for h in lags
    ]
end

lags = 1:30

ext_obs = [extremogram(data, u, lags) for u in thresholds]                      # one vector per threshold
ext_rep = [
    reduce(hcat, [extremogram(view(xrep, :, j), u, lags) for j in 1:M])         # lags x M
    for u in thresholds
]

ext_panels = map(eachindex(thresholds)) do k
    p = band_plot(collect(lags), ext_obs[k], ext_rep[k];
        xlabel = "Lag h", ylabel = "Extremogram",
        title = "Threshold q$(threshold_probs[k])",
        legend = (k == 1 ? :topright : false))
    # level expected under independence: P(X > u) = 1 - p
    hline!(p, [1 - threshold_probs[k]]; linestyle = :dot, color = :gray,
           label = (k == 1 ? "Independence level" : ""))
    p
end
p_extremogram = plot(ext_panels...; layout = (length(thresholds), 1),
                     size = (800, 300 * length(thresholds)))
display(p_extremogram)

