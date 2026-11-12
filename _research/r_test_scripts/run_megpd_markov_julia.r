library(JuliaCall)
library(tidyverse)

julia_setup()

julia_source("megpd/code/src/MMEGPD.jl")

n <- 10000
kappa_val <- 2
sigma_val <- 1
xi_val <- 0.5
delta <- julia_eval(":delta08")

sim <- julia_call("MMEGPD.simulate_megpd_chain", as.integer(n), kappa_val, sigma_val, xi_val, delta)

x <- sim

print(acf(x, lag.max = 100))
print(pacf(x, lag.max = 100))

library(posterior)
ess_bulk(x)
ess_tail(x)

library(egpd)
fitegpd(x, type = 1, family = "egpd")

# Quantiles simulation
julia_source("megpd/code/quantile_megpd_markov.jl")

m <- 50
quants <- as.vector(seq(0.001, 0.999, length.out = 400))

sim_qq <- julia_call("simulate_megpd_quantiles",
                        as.integer(m),
                        as.integer(n),
                        quants,
                        kappa_val,
                        sigma_val,
                        xi_val,
                        delta
)

Q <- sim_qq$Q

p <- quants

# Q <- read.csv("megpd/res/sim_qq_M50_n10k_xi05.csv")

theoretical <- egpd::qegpd(
    p,
    kappa = kappa_val,
    sigma = sigma_val,
    xi = xi_val
)

df <- tibble(
    p = p,
    theo = theoretical,
    lower = apply(Q, 2, quantile, 0.025, na.rm = TRUE),
    mean = apply(Q, 2, mean, na.rm = TRUE),
    median = apply(Q, 2, quantile, 0.50, na.rm = TRUE),
    upper = apply(Q, 2, quantile, 0.975, na.rm = TRUE)
)

df_spaghetti <- tibble(
    replicate = rep(seq_len(nrow(Q)), each = ncol(Q)),
    theo       = rep(df$theo, times = nrow(Q)),
    empirical  = c(t(Q))
)
p_qq_spaghetti <- ggplot() +

    geom_line(
        data = df_spaghetti,
        aes(
            x = theo,
            y = empirical,
            group = replicate
        ),
        colour = "lightblue",
        alpha = 0.7,
        linewidth = 0.3
    ) +

    geom_line(
        data = df,
        aes(theo, mean),
        colour = "blue",
        linewidth = 1
    ) +

    geom_abline(
        slope = 1,
        intercept = 0,
        colour = "red",
        linetype = 2
    ) +

    scale_x_log10() +
    scale_y_log10() +

    labs(
        title = "Radius Quantile distribution",
        x = "Theoretical EGPD quantile",
        y = "Empirical quantile"
    ) +

    theme_bw()

ggsave("megpd/figures/julia_spaghetti_radius_quantiles_deltastrongupper_M50_n1k_xi05.png", p_qq_spaghetti)
