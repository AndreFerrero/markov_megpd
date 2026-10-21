library(JuliaCall)
library(tidyverse)
source("megpd/code/megpd_fit.r")

extract_radial_estimates <- function(sims){
    fits <- list()
    r_param <- matrix(NA, M, 3)
    for (m in seq_len(M)){
        fits[[m]] <- fit_megpd_ts(sims[,m])$radial_fit
        r_param[m,] <- fits[[m]]$estimate
    }
    colnames(r_param) <- names(fits[[1]]$estimate)

    r_param
}

julia_setup()

julia_source("megpd/code/markov_megpd.jl")
julia_source("megpd/code/delta.jl")

n <- 10000
kappa_val <- 2
sigma_val <- 1
xi_val <- 0.5
delta <- julia_eval(":delta_strong_upper")

M <- 100

sims_delta_strong_upper <- julia_call("multiple_megpd_chains", as.integer(M), as.integer(n), kappa_val, sigma_val, xi_val, delta)

save(sims_delta_strong_upper, file = "megpd/res/sims_delta_strong_upper_n10k_M100_kappa2sigma1xi05.Rdata")
# load("megpd/res/sims_delta_strong_upper_n1k_M500_kappa2sigma1xi05.Rdata")

r_param_delta_strong_upper <- extract_radial_estimates(sims_delta_strong_upper)

round(colMeans(r_param_delta_strong_upper), 2)

source("code/models/copula_markov/copula_markov_model.r")
source("code/models/margins/egpd.r")
source("code/models/copulas/gumbel.r")

egpd_gumbel_model <- make_copula_markov_model(margin_egpd, copula_gumbel, stan_mod = NULL)

set.seed(123)
margin_param <- c(kappa = kappa_val, sigma = sigma_val, xi = xi_val)
copula_param <- 3

gumbel_sims <- matrix(NA, ncol = M, nrow = n)

for (m in seq_len(M)) {
    print(paste("Simulation: ", m))
    gumbel_sims[,m] <- egpd_gumbel_model$simulate(
        n = n,
        margin_param = margin_param,
        copula_param = copula_param
    )$x
}
save(gumbel_sims, file = "megpd/res/gumbel_sims_n1k_M500_kappa2sigma1xi05.Rdata")


r_param_gumbel <- extract_radial_estimates(gumbel_sims)

round(colMeans(r_param_gumbel), 3)
round(apply(r_param_gumbel, 2, "sd"), 3)

boxplot(
  r_param_gumbel,
  names = names(true_vals),
  col = c("skyblue", "lightgreen", "salmon"),
  ylab = "Estimate"
)

points(1:3, c(sigma_val, xi_val, kappa_val), pch = 19, col = "red")
