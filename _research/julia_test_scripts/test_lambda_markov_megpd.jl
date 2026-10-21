
using MMEGPD

using Random
using BSplines
using Optim
using Plots
using JLD2

Random.seed!(1234)

# ============================================================
# 1. Simulate data
# ============================================================

n = 5000

data = MMEGPD.simulate_megpd_chain(
    n,
    2,
    1,
    0.2,
    :delta_strong_upper
)

x = data


# ============================================================
# 2. Construct radial observations and spline basis
# ============================================================

R = x[1:end-1] .+ x[2:end]

r_grid = range(
    minimum(R),
    maximum(R),
    length=500
)

K = 10

basis, X =
    MMEGPD.bspline_basis(
        R;
        nbasis=K,
        degree=3
    )

B_grid =
    basismatrix(
        basis,
        collect(r_grid)
    )


# ============================================================
# 3. True delta function
# ============================================================

delta_true =
    [MMEGPD.delta_strong_upper(r) for r in r_grid]


# ============================================================
# 4. Lambda values to compare
# ============================================================

lambda_grid = [
    10.0,
    15.0,
    20.0,
    30.0,
    50.0,
    70.0,
    100.0,
    150.0,
    200.0,
    300.0,
    500.0,
    1000.0
]

println("Lambda grid:")
println(lambda_grid)


# ============================================================
# 5. Fit model for every lambda
# ============================================================

fits = Dict{Float64, Any}()
delta_hats = Dict{Float64, Vector{Float64}}()
objective_values = Dict{Float64, Float64}()

theta_hats = Dict{Float64, Vector{Float64}}()
beta_hats = Dict{Float64, Vector{Float64}}()

for lambda in lambda_grid

    println()
    println("==========================================")
    println("Fitting lambda = ", lambda)
    println("==========================================")

    fit =
        MMEGPD.fit_markov_megpd_splines(
            x,
            K;
            p = 2,
            lambda=lambda
        )

    fits[lambda] = fit

    # --------------------------------------------------------
    # Store parameters
    # --------------------------------------------------------

    theta_hats[lambda] = collect(fit.theta)
    beta_hats[lambda] = collect(fit.beta)

    # --------------------------------------------------------
    # Fitted delta on grid
    # --------------------------------------------------------

    log_delta_hat =
        B_grid * fit.beta

    delta_hat =
        exp.(log_delta_hat)

    delta_hats[lambda] =
        delta_hat

    # --------------------------------------------------------
    # Objective
    # --------------------------------------------------------

    objective_values[lambda] =
        fit.result.minimum

    println("theta = ", fit.theta)
    println("beta  = ", fit.beta)
    println("objective = ", fit.result.minimum)

end

save_file = "megpd/res/lambda_loop_n5k_p2.jld2"
@save save_file lambda_grid fits theta_hats beta_hats delta_hats objective_values r_grid delta_true x K n

println("Results saved to: ", save_file)

data = load(save_file)

lambda_grid = data["lambda_grid"]
fits = data["fits"]
theta_hats = data["theta_hats"]
beta_hats = data["beta_hats"]
delta_hats = data["delta_hats"]
objective_values = data["objective_values"]
r_grid = data["r_grid"]
delta_true = data["delta_true"]
x = data["x"]
K = data["K"]
n = data["n"]

# ============================================================
# 6. Plot fitted delta curves
# ============================================================

p1 = plot(
    r_grid,
    delta_true,
    linewidth=4,
    label="true δ(r)",
    xlabel="r",
    ylabel="δ(r)",
    title="Effect of smoothing parameter",
    ylims=(0, 1)
)

for lambda in lambda_grid

    plot!(
        p1,
        r_grid,
        delta_hats[lambda],
        linewidth=1.5,
        label="λ = $(round(lambda, sigdigits=3))"
    )

end

display(p1) 

kappa_hat = [theta_hats[λ][1] for λ in lambda_grid]
sigma_hat = [theta_hats[λ][2] for λ in lambda_grid]
xi_hat    = [theta_hats[λ][3] for λ in lambda_grid]

p_kappa = plot(
    lambda_grid,
    kappa_hat,
    xscale=:log10,
    marker=:circle,
    linewidth=2,
    xlabel="λ",
    ylabel="κ̂",
    title="Effect of penalisation on κ"
)

p_sigma = plot(
    lambda_grid,
    sigma_hat,
    xscale=:log10,
    marker=:circle,
    linewidth=2,
    xlabel="λ",
    ylabel="σ̂",
    title="Effect of penalisation on σ"
)

p_xi = plot(
    lambda_grid,
    xi_hat,
    xscale=:log10,
    marker=:circle,
    linewidth=2,
    xlabel="λ",
    ylabel="ξ̂",
    title="Effect of penalisation on ξ"
)

display(p_kappa)
display(p_sigma)
display(p_xi)