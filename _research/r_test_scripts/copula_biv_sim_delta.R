library(copula)
library(egpd)

source("C:/Users/andre/markov_megpd/_research/r_test_scripts/megpd_2stepfit.r")
source("C:/Users/andre/markov_megpd/_research/r_test_scripts/delta_copula.r")


# ============================================================
# 1. SIMULATION SETTINGS
# ============================================================

set.seed(123)

n <- 2000
B <- 200

xi_true <- 0.2
kappa_true <- 2
sigma_true <- 1

theta <- 4
cop <- gumbelCopula(theta)

# Output directory
output_dir <- "C:/Users/andre/markov_megpd/_research/data"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. HELPER FUNCTION: NUMERICAL INTEGRATION
# ============================================================

# Trapezoidal integration.
# x must contain distinct, finite values.

trapz <- function(x, y) {

  keep <- is.finite(x) & is.finite(y)

  x <- x[keep]
  y <- y[keep]

  if (length(x) < 2) {
    return(NA_real_)
  }

  ord <- order(x)

  x <- x[ord]
  y <- y[ord]

  # Remove duplicated x values, if any
  keep <- !duplicated(x)

  x <- x[keep]
  y <- y[keep]

  if (length(x) < 2) {
    return(NA_real_)
  }

  sum(
    diff(x) *
      (head(y, -1) + tail(y, -1)) / 2
  )
}


# ============================================================
# 3. STORAGE
# ============================================================

# Keep the entire fitted object for every repetition
fits <- vector("list", B)

# Keep the simulated observations
datasets <- vector("list", B)

# Keep estimated and theoretical delta curves
delta_curves <- vector("list", B)

# Track fitting errors
fit_errors <- rep(NA_character_, B)
delta_errors <- rep(NA_character_, B)

# Numerical results
metrics <- data.frame(
  iteration = seq_len(B),

  xi_hat = rep(NA_real_, B),

  radial_min = rep(NA_real_, B),
  radial_max = rep(NA_real_, B),
  radial_width = rep(NA_real_, B),

  ISE = rep(NA_real_, B),
  IRMSE = rep(NA_real_, B),
  relative_ISE = rep(NA_real_, B),
  relative_IRMSE = rep(NA_real_, B),

  n_curve_points = rep(NA_integer_, B),

  fit_success = rep(FALSE, B),
  delta_success = rep(FALSE, B)
)


# ============================================================
# 4. SIMULATION LOOP
# ============================================================

for (b in seq_len(B)) {

  message("Simulation ", b, " / ", B)

  # ----------------------------------------------------------
  # Generate a dataset
  # ----------------------------------------------------------

  U <- rCopula(n, cop)

  X <- apply(
    U,
    2,
    qegpd,
    kappa = kappa_true,
    sigma = sigma_true,
    xi = xi_true
  )

  datasets[[b]] <- X

  R <- X[, 1] + X[, 2]


  # ----------------------------------------------------------
  # Fit MEGPD
  # ----------------------------------------------------------

  fit <- tryCatch(

    fit_megpd_bivariate(
      X[, 1],
      X[, 2]
    ),

    error = function(e) {

      fit_errors[b] <<- conditionMessage(e)

      message(
        "Fit failed at iteration ", b, ": ",
        conditionMessage(e)
      )

      NULL
    }
  )

  if (is.null(fit)) {
    next
  }

  # Save the complete fitted object
  fits[[b]] <- fit

  metrics$fit_success[b] <- TRUE


  # ----------------------------------------------------------
  # Extract the radial shape estimate
  # ----------------------------------------------------------

  radial_est <- fit$radial_fit$estimate

  if (!is.null(radial_est) &&
      "xi" %in% names(radial_est)) {

    metrics$xi_hat[b] <- unname(radial_est["xi"])

  }


  # ----------------------------------------------------------
  # Estimate delta on the dataset-specific radial domain
  # ----------------------------------------------------------

  delta_est <- tryCatch(

    plot_delta(
      fit$angular_fit,
      R
    ),

    error = function(e) {

      delta_errors[b] <<- conditionMessage(e)

      message(
        "Delta estimation failed at iteration ", b, ": ",
        conditionMessage(e)
      )

      NULL
    }
  )

  if (is.null(delta_est)) {
    next
  }


  # ----------------------------------------------------------
  # Calculate the theoretical delta on exactly the same grid
  # ----------------------------------------------------------

  delta_true <- tryCatch(

    delta_copula(
      r = delta_est$R,
      cop = cop,
      kappa = kappa_true,
      sigma = sigma_true,
      xi = xi_true
    ),

    error = function(e) {

      delta_errors[b] <<- conditionMessage(e)

      message(
        "Theoretical delta failed at iteration ", b, ": ",
        conditionMessage(e)
      )

      NULL
    }
  )

  if (is.null(delta2_true)) {
    next
  }

  # ----------------------------------------------------------
  # Store all curve information
  # ----------------------------------------------------------

  curve <- data.frame(
    R = delta_est$R,
    delta_hat = delta_est$delta_hat,
    delta_lo = delta_est$delta_lo,
    delta_hi = delta_est$delta_hi,
    delta_true = as.numeric(delta_true)
  )

  # Remove invalid rows and sort by radius
  curve <- curve[
    is.finite(curve$R) &
      is.finite(curve$delta_hat) &
      is.finite(curve$delta_true),
  ]

  curve <- curve[order(curve$R), ]

  # Remove duplicated radial values
  curve <- curve[!duplicated(curve$R), ]

  delta_curves[[b]] <- curve


  # ----------------------------------------------------------
  # Calculate numerical accuracy metrics
  # ----------------------------------------------------------

  if (nrow(curve) < 2) {
    next
  }

  r <- curve$R

  delta_hat <- curve$delta_hat
  delta_true <- curve$delta_true

  error_squared <- (delta_hat - delta_true)^2
  true_squared <- delta_true^2

  # Dataset-specific radial domain
  r_min <- min(r)
  r_max <- max(r)
  r_width <- r_max - r_min

  if (!is.finite(r_width) || r_width <= 0) {
    next
  }

  # Integrated squared error
  ISE <- trapz(
    r,
    error_squared
  )

  # Integrated mean squared error
  IMSE <- ISE / r_width

  # Integrated root mean squared error
  IRMSE <- sqrt(IMSE)

  # Relative integrated squared error
  denominator <- trapz(
    r,
    true_squared
  )

  relative_ISE <- if (
    is.finite(denominator) && denominator > 0
  ) {
    ISE / denominator
  } else {
    NA_real_
  }

  # Relative integrated root mean squared error
  relative_IRMSE <- sqrt(relative_ISE)


  # ----------------------------------------------------------
  # Save metrics for this repetition
  # ----------------------------------------------------------

  metrics$radial_min[b] <- r_min
  metrics$radial_max[b] <- r_max
  metrics$radial_width[b] <- r_width

  metrics$ISE[b] <- ISE
  metrics$IRMSE[b] <- IRMSE

  metrics$relative_ISE[b] <- relative_ISE
  metrics$relative_IRMSE[b] <- relative_IRMSE

  metrics$n_curve_points[b] <- nrow(curve)

  metrics$delta_success[b] <- TRUE


  # ----------------------------------------------------------
  # Periodic checkpoint
  # ----------------------------------------------------------

  if (b %% 25 == 0) {

    saveRDS(
      list(
        settings = list(
          n = n,
          B = B,
          xi_true = xi_true,
          kappa_true = kappa_true,
          sigma_true = sigma_true,
          copula = cop
        ),

        fits = fits,
        datasets = datasets,
        delta_curves = delta_curves,
        metrics = metrics,
        fit_errors = fit_errors,
        delta_errors = delta_errors
      ),

      file = file.path(
        output_dir,
        "megpd_simulation_checkpoint.rds"
      )
    )

    write.csv(
      metrics,
      file.path(output_dir, "delta_metrics_checkpoint.csv"),
      row.names = FALSE
    )
  }
}


# ============================================================
# 5. SAVE FINAL RESULTS
# ============================================================

simulation_results <- list(

  settings = list(
    n = n,
    B = B,
    xi_true = xi_true,
    kappa_true = kappa_true,
    sigma_true = sigma_true,
    copula = cop,
    theta = theta
  ),

  fits = fits,

  datasets = datasets,

  delta_curves = delta_curves,

  metrics = metrics,

  fit_errors = fit_errors,

  delta_errors = delta_errors
)


saveRDS(
  simulation_results,
  file = file.path(
    output_dir,
    "copula_biv_sim.rds"
  )
)


write.csv(
  metrics,
  file = file.path(
    output_dir,
    "copula_biv_delta.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 6. SUMMARY OF RADIAL XI ESTIMATES
# ============================================================

valid_xi <- metrics$xi_hat[
  is.finite(metrics$xi_hat)
]

summary_xi <- c(

  repetitions = B,

  successful_fits = sum(metrics$fit_success),

  successful_delta = sum(metrics$delta_success),

  mean_xi = if (length(valid_xi) > 0) {
    mean(valid_xi)
  } else {
    NA_real_
  },

  sd_xi = if (length(valid_xi) > 1) {
    sd(valid_xi)
  } else {
    NA_real_
  },

  median_xi = if (length(valid_xi) > 0) {
    median(valid_xi)
  } else {
    NA_real_
  }
)

print(summary_xi)


# ============================================================
# 7. SUMMARY OF DELTA ACCURACY
# ============================================================

valid_metrics <- metrics[
  metrics$delta_success,
]

if (nrow(valid_metrics) > 0) {

  summary_delta <- data.frame(

    metric = c(
      "ISE",
      "IRMSE",
      "relative_ISE",
      "relative_IRMSE"
    ),

    mean = c(
      mean(valid_metrics$ISE, na.rm = TRUE),
      mean(valid_metrics$IRMSE, na.rm = TRUE),
      mean(valid_metrics$relative_ISE, na.rm = TRUE),
      mean(valid_metrics$relative_IRMSE, na.rm = TRUE)
    ),

    median = c(
      median(valid_metrics$ISE, na.rm = TRUE),
      median(valid_metrics$IRMSE, na.rm = TRUE),
      median(valid_metrics$relative_ISE, na.rm = TRUE),
      median(valid_metrics$relative_IRMSE, na.rm = TRUE)
    ),

    sd = c(
      sd(valid_metrics$ISE, na.rm = TRUE),
      sd(valid_metrics$IRMSE, na.rm = TRUE),
      sd(valid_metrics$relative_ISE, na.rm = TRUE),
      sd(valid_metrics$relative_IRMSE, na.rm = TRUE)
    )
  )

  print(summary_delta)

  write.csv(
    summary_delta,
    file.path(output_dir, "copula_biv_delta_summary"),
    row.names = FALSE
  )
}


# ============================================================
# 8. DISTRIBUTION OF XI ESTIMATES
# ============================================================

if (length(valid_xi) > 0) {

  hist(
    valid_xi,
    breaks = 30,
    main = expression("Distribution of estimated " * hat(xi)),
    xlab = expression(hat(xi)),
    col = "lightblue",
    border = "white"
  )

  abline(
    v = mean(valid_xi),
    col = "red",
    lwd = 2
  )

  legend(
    "topright",
    legend = "Mean estimate",
    col = "red",
    lty = 1,
    lwd = 2,
    bty = "n"
  )
}


# ============================================================
# 9. PLOT ALL ESTIMATED AND THEORETICAL DELTA CURVES
# ============================================================

pdf(
  file.path(output_dir, "delta_curve_comparisons.pdf"),
  width = 8,
  height = 6
)

for (b in seq_len(B)) {

  curve <- delta_curves[[b]]

  if (is.null(curve) || nrow(curve) < 2) {
    next
  }

  y_range <- range(
    c(
      0,
      curve$delta_hat,
      curve$delta_lo,
      curve$delta_hi,
      curve$delta_true
    ),
    finite = TRUE
  )

  plot(
    curve$R,
    curve$delta_hat,
    type = "l",
    lwd = 2,
    col = "black",
    xlab = "R",
    ylab = expression(delta(R)),
    ylim = y_range,
    main = paste("Simulation", b)
  )

  lines(
    curve$R,
    curve$delta_lo,
    col = "red",
    lty = 2
  )

  lines(
    curve$R,
    curve$delta_hi,
    col = "red",
    lty = 2
  )

  lines(
    curve$R,
    curve$delta_true,
    col = "blue",
    lwd = 2
  )

  legend(
    "topright",
    legend = c(
      "Estimated delta",
      "95% confidence interval",
      "Theoretical delta"
    ),
    col = c("black", "red", "blue"),
    lty = c(1, 2, 1),
    lwd = c(2, 1, 2),
    bty = "n"
  )
}

dev.off()


# ============================================================
# 10. PLOT DISTRIBUTION OF DELTA ERRORS
# ============================================================

if (nrow(valid_metrics) > 0) {

  par(mfrow = c(1, 2))

  hist(
    valid_metrics$IRMSE,
    breaks = 25,
    main = "Integrated RMSE",
    xlab = "IRMSE",
    col = "lightblue",
    border = "white"
  )

  hist(
    valid_metrics$relative_IRMSE,
    breaks = 25,
    main = "Relative integrated RMSE",
    xlab = "Relative IRMSE",
    col = "lightblue",
    border = "white"
  )

  par(mfrow = c(1, 1))
}


# ============================================================
# 11. COMPLETION
# ============================================================

message("Simulation completed.")

message(
  "Successful fits: ",
  sum(metrics$fit_success),
  " / ",
  B
)

message(
  "Successful delta comparisons: ",
  sum(metrics$delta_success),
  " / ",
  B
)

message("Results saved to: ", output_dir)