library(JuliaCall)
library(tidyverse)
library(posterior)
library(egpd)

julia_setup()

julia_source("megpd/code/src/MMEGPD.jl")


# ============================================================
# PARAMETERS
# ============================================================

n <- 10000
n_chains <- 100

kappa_val <- 2
sigma_val <- 1
xi_val <- 0.5

burn_in_prop <- 0.2

delta <- julia_eval(":delta08")


# ============================================================
# SIMULATE MULTIPLE CHAINS
# ============================================================

chains <- julia_call(
    "MMEGPD.multiple_megpd_chains",
    as.integer(n_chains),
    as.integer(n),
    kappa_val,
    sigma_val,
    xi_val,
    delta,
    x0 = 1.0,
    burn_in_prop = burn_in_prop
)

dim(chains)

# ============================================================
# FIT EGPD TO EACH CHAIN
# ============================================================

estimates <- matrix(
    NA_real_,
    nrow = n_chains,
    ncol = 3
)

colnames(estimates) <- c(
    "sigma",
    "xi",
    "kappa"
)

for (j in seq_len(n_chains)) {

    cat("Fitting chain", j, "of", n_chains, "\n")

    x <- chains[, j]

    fit <- tryCatch(
        egpd::fitegpd(
            x,
            type = 1,
            family = "egpd"
        ),
        error = function(e) NULL
    )

    if (!is.null(fit)) {

        # Replace these three lines with the appropriate
        # extraction from your fitegpd() object
        estimates[j, "sigma"] <- fit$estimate[1]
        estimates[j, "xi"] <- fit$estimate[2]
        estimates[j, "kappa"]    <- fit$estimate[3]
    }
}

colMeans(estimates, na.rm = TRUE)

apply(
    estimates,
    2,
    sd,
    na.rm = TRUE
)