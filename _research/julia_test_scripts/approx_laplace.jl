
using MMEGPD
using Random
using DataFrames
using Base.Threads

Random.seed!(1234)


x = MMEGPD.simulate_megpd_chain(2000, 2, 1, 0.2, :delta_strong_upper)

R = x[1:(end-1)] .+ x[2:end]

r_grid = range(
    minimum(R),
    maximum(R),
    length=500
)

K = 3


# lambdas = 10 .^ [1, 2, 3]


# # ---------------------------------------------------------
# # 3. Fit model for each lambda
# # ---------------------------------------------------------

# laplace_values = Float64[]
# numerical_values = Float64[]
# numerical_errors = Float64[]

# fits = []

# for lambda in lambdas

#     println("lambda = ", lambda)

#     fit = MMEGPD.fit_markov_megpd_splines(
#         x,
#         K;
#         lambda=lambda,
#         spline=:crspline,
#         knot_method=:quantile,
#         map=true,
#         verbose=false
#     )

#     push!(fits, fit)

#     # Laplace approximation
#     L_lap = MMEGPD.laplace_log_marginal(
#         fit;
#         drop_constant=false
#     )

#     println("Starting Integration for lambda = $lambda")
#     # Numerical integration
#     L_num, err = MMEGPD.numerical_log_marginal(
#         fit;
#         x=x,
#         se_w=5.0,
#         drop_constant=false,
#         rtol=1e-1,
#         atol=1e-5
#     )

#     push!(laplace_values, L_lap)
#     push!(numerical_values, L_num)
#     push!(numerical_errors, err)

# end

# results = DataFrame(
#     lambda=lambdas,
#     laplace=laplace_values,
#     numerical=numerical_values,
#     difference=laplace_values .- numerical_values,
#     numerical_error=numerical_errors
# )

# println(results)

# Increasing withs
w = [1.0, 2.0, 3.0, 4.0, 5.0]
lambda = 10.0

fit = MMEGPD.fit_markov_megpd_splines(
        x,
        K;
        lambda=lambda,
        spline=:crspline,
        knot_method=:quantile,
        map=true,
        verbose=false
    )

# Laplace approximation
L_lap = MMEGPD.laplace_log_marginal(
    fit;
    drop_constant=false
)

n = length(w)

numerical_values = Vector{Float64}(undef, n)
numerical_errors = Vector{Float64}(undef, n)

println("Starting loop")

Threads.@threads for i in eachindex(w)

    se_w = w[i]

    Core.println("Starting integration with se_w = $se_w")

    L_num, err = MMEGPD.numerical_log_marginal(
        fit;
        x=x,
        se_w=se_w,
        drop_constant=false,
        rtol=1e-1,
        atol=1e-5
    )

    numerical_values[i] = L_num
    numerical_errors[i] = err

    Core.println("Finished integration with se_w = $se_w")
end

results = DataFrame(
    se=w,
    laplace=fill(L_lap, n),
    numerical=numerical_values,
    difference=L_lap .- numerical_values,
    numerical_error=numerical_errors
)

println(results)