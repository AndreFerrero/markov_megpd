
using MMEGPD

using Random
using BSplines
using Optim
using Plots
using JLD2

Random.seed!(1234)

n = 10000

data = MMEGPD.simulate_megpd_chain(
    n,
    2,
    1,
    0.2,
    :delta_strong_upper
)

N_train = floor(Int, 0.5 * length(data))

x_train = data[1:N_train]
x_test  = data[N_train+1:end]

lambda_grid = [
    10.0,
    50.0,
    100.0,
    120.0,
    150.0,
    200.0,
    300.0,
    500.0,
    1000.0
]

K = 10
res_eval = MMEGPD.select_lambda(
    x_train,
    x_test,
    K,
    lambda_grid;
    p = 2
)

using JLD2

@save "megpd/res/lambda_eval_K20.jld2" res_eval