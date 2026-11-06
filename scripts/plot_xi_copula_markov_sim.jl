using DrWatson
@quickactivate "markov_megpd"

using DataFrames
using Statistics
using StatsPlots
using Random

include(srcdir("copulasim_funs.jl"))

name = isempty(ARGS) ? first(readdir(datadir(FITS_SUBDIR))) : ARGS[1]

df = collect_results(datadir("copulasimfits/$name"))

# Select relevant columns and remove missing values
xi_results = dropmissing(
    select(df, :family, :cop_param, :true_xi, :xi_hat, :rep)
)

# ---------------------------------------------------------------------------
# Plot settings
# ---------------------------------------------------------------------------

families = ["gumbel", "joe", "normal"]

cop_params = Dict(
    "gumbel" => [2.0, 4.0],
    "joe"    => [2.0, 4.0],
    "normal" => [0.5, 0.75],
)

family_colors = Dict(
    "gumbel" => :steelblue,
    "joe"    => :orange,
    "normal" => :seagreen,
)

# ---------------------------------------------------------------------------
# Function to create a boxplot for a given true xi
# ---------------------------------------------------------------------------

function plot_xi(data, xi_value)

    plot_data = filter(
        row -> isapprox(row.true_xi, xi_value; atol=1e-8),
        data
    )

    p = plot(
        xlabel = "Copula scenario",
        ylabel = "Estimated ξ",
        title = "Estimation of ξ: true ξ = $xi_value",
        legend = :topright,
        xticks = (1:6, String[]),
        xtickfont = font(10),
        guidefont = font(12),
        titlefont = font(13),
        size = (950, 600),
        framestyle = :box,
        grid = :y,
    )

    labels = String[]
    position = 0

    Random.seed!(1234)

    for family in families
        for cop_param in cop_params[family]

            position += 1

            scenario_data = filter(
                row -> row.family == family &&
                       isapprox(row.cop_param, cop_param; atol=1e-8),
                plot_data
            )

            estimates = scenario_data.xi_hat

            push!(
                labels,
                family == "normal" ?
                    "Gaussian\nρ = $cop_param" :
                    "$(uppercasefirst(family))\nθ = $cop_param"
            )

            # Boxplot
            boxplot!(
                p,
                fill(position, length(estimates)),
                estimates,
                color = family_colors[family],
                alpha = 0.65,
                linewidth = 1.5,
                outlier = false,
                label = "",
                bar_width = 0.6,
            )

            # Individual Monte Carlo estimates
            jitter = rand(length(estimates)) .* 0.18 .- 0.09

            scatter!(
                p,
                position .+ jitter,
                estimates,
                color = family_colors[family],
                markersize = 3,
                alpha = 0.65,
                markerstrokewidth = 0,
                label = "",
            )
        end
    end

    # True xi reference line
    hline!(
        p,
        [xi_value],
        color = :red,
        linestyle = :dash,
        linewidth = 2,
        label = "True ξ = $xi_value",
    )

    # Scenario labels
    plot!(
        p,
        xticks = (1:6, labels),
    )

    return p
end

# ---------------------------------------------------------------------------
# Generate the two plots
# ---------------------------------------------------------------------------

p_xi_neg = plot_xi(xi_results, -0.2)
p_xi_pos = plot_xi(xi_results, 0.2)

display(p_xi_neg)
display(p_xi_pos)

# ---------------------------------------------------------------------------
# Save plots
# ---------------------------------------------------------------------------

# savefig(p_xi_neg, plotsdir("xi_boxplot_minus_02.png"))
# savefig(p_xi_pos, plotsdir("xi_boxplot_plus_02.png"))