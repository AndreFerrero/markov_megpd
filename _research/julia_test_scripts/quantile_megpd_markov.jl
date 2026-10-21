using Statistics
using StatsBase
using Base.Threads
using ProgressMeter
using CSV
using Tables
using Random

using MMEGPD

"""
Run M replications of the MEGPD chain simulation and compute quantiles of Rt = x[t] + x[t-1].
"""
function simulate_megpd_quantiles(
    M,
    n_steps,
    probs,
    kappa,
    sigma,
    xi,
    delta_name::Symbol;
    x0=1.0,
    burn_in_prop=0.0,
)
    Random.seed!(1234)  # Set a seed for reproducibility

    P = length(probs)

    Q = Matrix{Float64}(undef, M, P)
    success = falses(M)
    error_message = Vector{Union{String,Nothing}}(undef, M)

    for m in 1:M
        println("\n Replication $m of $M...")
        
        try
            sim = MMEGPD.simulate_megpd_chain(
                n_steps,
                kappa,
                sigma,
                xi,
                delta_name;
                x0 = x0,
                burn_in_prop = burn_in_prop,
                show_progress = false
            )

            x = sim

            # Rt = x_t + x_{t-1}
            Rt = x[2:end] .+ x[1:(end-1)]

            # quantiles
            Q[m, :] .= quantile(Rt, probs)

            success[m] = true
            error_message[m] = nothing

        catch e
            success[m] = false
            error_message[m] = sprint(showerror, e)
            @info "FAILED replication $m" exception=e
        end
    end

    return (
        probs=probs,
        Q=Q,
        success=success,
        error_message=error_message
    )
end