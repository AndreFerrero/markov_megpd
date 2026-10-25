using DrWatson
@quickactivate "markov_megpd"

using DataFrames
using Statistics
using RCall

include(srcdir("copulasim_funs.jl"))

name = isempty(ARGS) ? first(readdir(datadir(FITS_SUBDIR))) : ARGS[1]

# df = collect_results(datadir("copulasim"))
load_simulated_data(name)

df = collect_results(datadir("copulasimfits/$name"))

# Select relevant columns and remove missing values
xi_results = dropmissing(
    select(df, :family, :tau, :true_xi, :xi_hat, :rep)
)

# Pass the Julia DataFrame to R
@rput xi_results

R"""
library(ggplot2)
library(dplyr)

# Prepare the data
xi_results <- xi_results %>%
    mutate(
        family = factor(
            family,
            levels = c("gumbel", "joe", "normal")
        ),
        tau = factor(
            tau,
            levels = c(0.5, 0.75),
            labels = c("tau = 0.5", "tau = 0.75")
        )
    )

# Function to create a boxplot for a given true xi
plot_xi <- function(data, xi_value) {

    plot_data <- data %>%
        filter(true_xi == xi_value)

    ggplot(
        plot_data,
        aes(
            x = interaction(family, tau, sep = "\n"),
            y = xi_hat,
            fill = family
        )
    ) +
        geom_boxplot(
            width = 0.6,
            alpha = 0.7,
            outlier.shape = NA
        ) +
        geom_jitter(
            aes(color = family),
            width = 0.12,
            size = 2,
            alpha = 0.8
        ) +
        geom_hline(
            yintercept = xi_value,
            color = "red",
            linetype = "dashed",
            linewidth = 0.9
        ) +
        scale_fill_manual(
            values = c(
                "gumbel" = "#377EB8",
                "joe" = "#E69F00",
                "normal" = "#009E73"
            )
        ) +
        scale_color_manual(
            values = c(
                "gumbel" = "#377EB8",
                "joe" = "#E69F00",
                "normal" = "#009E73"
            ),
            guide = "none"
        ) +
        labs(
            title = paste("Estimation of xi: true xi =", xi_value),
            x = "Copula family and Kendall's tau",
            y = expression(hat(xi)),
            fill = "Copula family"
        ) +
        theme_bw(base_size = 13) +
        theme(
            plot.title = element_text(
                face = "bold",
                hjust = 0.5
            ),
            axis.text.x = element_text(
                angle = 0,
                hjust = 0.5
            ),
            legend.position = "bottom",
            panel.grid.minor = element_blank()
        )
}

# Plot for true xi = -0.2
p_xi_neg <- plot_xi(xi_results, -0.2)

# Plot for true xi = 0.2
p_xi_pos <- plot_xi(xi_results, 0.2)

# Display both plots
print(p_xi_neg)
print(p_xi_pos)

# Save plots
ggsave(
    "xi_boxplot_minus_02.png",
    plot = p_xi_neg,
    width = 10,
    height = 6,
    dpi = 300
)

ggsave(
    "xi_boxplot_plus_02.png",
    plot = p_xi_pos,
    width = 10,
    height = 6,
    dpi = 300
)
"""