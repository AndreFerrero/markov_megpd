"""
    intervals_estimator(T) -> θ̂

Intervals estimator of the extremal index (Ferro & Segers, 2003, eq. 4),
computed from the vector `T` of interexceedance times (T_i = S_{i+1} - S_i).
"""
function intervals_estimator(T::AbstractVector{<:Real})
    n = length(T)                       # n = N - 1
    n == 0 && return 1.0
    if maximum(T) <= 2
        θ = 2 * sum(T)^2 / (n * sum(abs2, T))
    else
        θ = 2 * sum(t -> t - 1, T)^2 / (n * sum(t -> (t - 1) * (t - 2), T))
    end
    return min(1.0, θ)
end

"""
    decluster_intervals(x, u; theta = nothing)

Interexceedance-times declustering (Ferro & Segers, 2003), automated version.

Arguments
- `x`     : observed series (vector, indexed by time 1:length(x))
- `u`     : threshold; exceedances are the times `i` with `x[i] > u`
- `theta` : extremal index to use. If `nothing` (default), it is estimated with
            the intervals estimator.

Method
1. Exceedance times S_1 < ... < S_N and interexceedance times T_i = S_{i+1} - S_i.
2. θ̂ from the intervals estimator (unless supplied).
3. Number of clusters N_C = 1 + floor(θ̂ (N - 1)) (capped at N).
4. The N_C - 1 largest interexceedance times are taken as the times separating
   clusters; every other gap stays inside a cluster. (Ties at the cut-off are
   broken by order of occurrence, so exactly N_C - 1 breaks are made.)

Returns a NamedTuple with
- `theta`      : extremal index used
- `nclusters`  : N_C
- `run_length` : smallest interexceedance time kept as a cluster separator
                 (equivalent run-length parameter; `Inf` if there is one cluster)
- `clusters`   : Vector of Vectors of indices (into `x`) of the exceedances per cluster
- `maxima_idx` : index in `x` of the maximum of each cluster
- `maxima`     : corresponding cluster maxima values
"""
function decluster_intervals(x::AbstractVector{<:Real}, u::Real; theta = nothing)
    S = findall(>(u), x)
    N = length(S)
    N == 0 && error("no exceedances over the threshold u = $u")

    if N == 1
        return (theta = 1.0, nclusters = 1, run_length = Inf,
                clusters = [S], maxima_idx = S, maxima = x[S])
    end

    T = diff(S)                                   # interexceedance times, length N-1
    θ = theta === nothing ? intervals_estimator(T) : float(theta)
    0 < θ <= 1 || error("theta must lie in (0, 1], got $θ")

    nc = min(N, 1 + floor(Int, θ * (N - 1)))      # number of clusters

    # positions k such that the gap between S[k] and S[k+1] separates two clusters
    order  = sortperm(T; rev = true)              # stable: ties by order of occurrence
    breaks = sort!(order[1:(nc - 1)])
    run_length = nc > 1 ? minimum(T[breaks]) : Inf

    clusters = Vector{Vector{Int}}()
    start = 1
    for b in breaks
        push!(clusters, S[start:b])
        start = b + 1
    end
    push!(clusters, S[start:end])

    maxima_idx = [idx[argmax(x[idx])] for idx in clusters]

    return (theta = θ, nclusters = nc, run_length = run_length,
            clusters = clusters, maxima_idx = maxima_idx, maxima = x[maxima_idx])
end
