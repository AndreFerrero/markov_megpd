using LinearAlgebra
using ForwardDiff

# ------------------------------------------------------------
# Compare two ways of computing log(det(H))
# ------------------------------------------------------------

using MMEGPD
using Random
using BSplines

Random.seed!(1234)


data = MMEGPD.simulate_megpd_chain(5000, 2, 1, 0.2, :delta_strong_upper)

x = data

R = x[1:end-1] .+ x[2:end]

r_grid =
    range(
        minimum(R),
        maximum(R),
        length=500
    )

K = 10
lam = 20

fit = MMEGPD.fit_markov_megpd_splines(x, K; lambda=lam, p = 2)


function compare_hessian_logdet(fit)

    params_hat = fit.params_hat
    negloglik  = fit.negloglik

    H = ForwardDiff.hessian(negloglik, params_hat)
    H = Symmetric((H + H') / 2)

    println("Hessian dimension: ", size(H))
    println()

    # --------------------------------------------------------
    # Method 1: Cholesky
    # --------------------------------------------------------

    println("=== Method 1: Cholesky ===")

    chol_success = false
    logdet_chol = NaN

    try
        C = cholesky(H)

        logdet_chol = 2 * sum(log.(diag(C.U)))
        chol_success = true

        println("Cholesky succeeded:      true")
        println("logdet(H) via Cholesky:  ", logdet_chol)

    catch e
        println("Cholesky succeeded:      false")
        println("Reason: ", e)
    end

    println()

    # --------------------------------------------------------
    # Method 2: Julia's logdet
    # --------------------------------------------------------

    println("=== Method 2: LinearAlgebra.logdet ===")

    logdet_builtin = logdet(H)

    println("logdet(H) via logdet:    ", logdet_builtin)

    println()

    # --------------------------------------------------------
    # Comparison
    # --------------------------------------------------------

    println("=== Comparison ===")

    if chol_success
        difference = logdet_chol - logdet_builtin

        println("Cholesky logdet:         ", logdet_chol)
        println("Built-in logdet:         ", logdet_builtin)
        println("Absolute difference:     ", abs(difference))
        println("Relative difference:     ",
                abs(difference) / max(abs(logdet_builtin), eps()))
    else
        println("Cannot compare: Cholesky failed.")
    end

    println()

    # --------------------------------------------------------
    # Extra diagnostic: eigenvalues
    # --------------------------------------------------------

    println("=== Hessian diagnostics ===")

    eigenvalues = eigvals(H)

    println("Minimum eigenvalue:      ", minimum(eigenvalues))
    println("Maximum eigenvalue:      ", maximum(eigenvalues))
    println("Condition number:        ", cond(H))

    println()
    println("Hessian positive definite: ",
            all(eigenvalues .> 0))

    return H
end

compare_hessian_logdet(fit)

using LinearAlgebra
using BenchmarkTools

# Hessian from your fitted model
H = ForwardDiff.hessian(fit.negloglik, fit.params_hat)
H = Symmetric((H + H') / 2)

println("Hessian dimension: ", size(H))

# Check that both methods agree before benchmarking
C = cholesky(H)

logdet_chol = 2 * sum(log.(diag(C.U)))
logdet_builtin = logdet(H)

println("\n=== Correctness check ===")
println("Cholesky:     ", logdet_chol)
println("LinearAlgebra: ", logdet_builtin)
println("Difference:   ", abs(logdet_chol - logdet_builtin))

println("\n=== Cholesky benchmark ===")
@btime begin
    C = cholesky($H)
    2 * sum(log.(diag(C.U)))
end

println("\n=== LinearAlgebra.logdet benchmark ===")
@btime logdet($H)