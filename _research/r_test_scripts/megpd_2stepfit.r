# Load required libraries
library(mgcv)
library(egpd) # Ensure your specific EGPD package is loaded

fit_megpd_bivariate <- function(x, y, bernstein = FALSE, m_degree = NULL, k = 10, type_s = "cr", mps = 2) {
  # 0. Data Preparation
  # The model is defined for positive intensities (e.g., river discharge)
  n <- length(x)
  R <- x + y # Radial component (sum-norm)
  V <- log(x / y) # Angular component (log-ratio)
  V <- V - mean(V) # Centering V to have mean 0 (as per the Z ~ zero-mean assumption)

  # 1. First Step: Radial component EGPD
  # The paper suggests m = floor(0.5 * n / log(n)) for Bernstein degree
  if (bernstein | !is.null(m_degree)) {
    if (is.null(m_degree)) {
      m_degree <- floor(0.5 * n / log(n))
    }

    fit_r <- egpd::fitegpd(
      R,
      type = 1,
      method = "bernstein",
      bernstein.m = m_degree
    )
  } else {
    fit_r <- egpd::fitegpd(R)
  }

  # 2. Second Step: Logistic-heteroscedastic modeling
  # We model V | R=r ~ N(0, delta(r)^2)
  # According to Eq (19), log(delta(r)) is a linear combination of basis functions.
  # We use mgcv's 'gauslss' family to model the scale (sigma) specifically.

  # In 'gauslss', the first formula is for the mean (mu), the second for log(sigma).
  # The paper assumes a zero-mean exchangeable vector Z (Eq 11).
  # Therefore, we fix intercept to 0 and no covariates for the mean.

  # K = 10 is the default basis dimension suggested in the paper (Section 4.1).
  fit_delta <- gam(
    list(
      V ~ 1, # Mean fixed at 0 (as per the Z ~ zero-mean assumption)
      ~ s(R, k = k, bs = type_s, m = mps) # log(sigma) modeled as a cubic spline of R
    ),
    family = gaulss(),
    data = data.frame(V = V, R = R)
  )

  # 3. Return Results
  results <- list(
    radial_fit = fit_r,
    angular_fit = fit_delta
  )

  return(results)
}

fit_megpd_ts <- function(x_ts, bernstein = FALSE, m_degree = NULL, k = 10, type_s = "cr", mps = 2) {
  x <- x_ts[-length(x_ts)]
  y <- x_ts[-1]

  fit_megpd_bivariate(x, y, bernstein, m_degree, k, type_s, mps)
}

plot_delta <- function(fit_delta, r_grid, plot = TRUE) {

  newdat <- data.frame(R = r_grid, V = 0)

  pred <- predict(fit_delta, newdata = newdat, type = "link", se.fit = TRUE)

  eta <- pred$fit[, 2]
  se <- pred$se.fit[, 2]
  eta_lo <- eta - 1.96 * se
  eta_hi <- eta + 1.96 * se

  sigma_hat <- exp(eta)
  sigma_lo <- exp(eta_lo)
  sigma_hi <- exp(eta_hi)

  if (plot) {
    plot(r_grid, sigma_hat,
      type = "l", lwd = 2,
      xlab = "R",
      ylab = expression(delta(R)),
      ylim = c(0, max(sigma_hi) * 1.2)
    )

    lines(r_grid, sigma_lo, col = "red", lty = 2)
    lines(r_grid, sigma_hi, col = "red", lty = 2)
  }

  invisible(data.frame(
    R = r_grid,
    delta_hat = sigma_hat,
    delta_lo = sigma_lo,
    delta_hi = sigma_hi
  ))
}
