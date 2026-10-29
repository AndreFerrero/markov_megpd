
library(copula)
library(egpd)

# Theoretical conditional variance:
# Var(log(X/Y) | X + Y = r)
#
# Assumes:
#   - Identical EGPD margins
#   - Exchangeable copula
#   - Radius R = X + Y
#   - Angular variable Z = log(X/Y)

delta_copula <- function(
    r,
    cop,
    kappa = 2,
    sigma = 1,
    xi = 0.2,
    rel.tol = 1e-7,
    subdivisions = 500L
) {

  # stopifnot(all(r > 0), length(r) >= 1)

  # EGPD marginal CDF and density
  Fm <- function(x) {
    pegpd(x, kappa = kappa, sigma = sigma, xi = xi)
  }

  fm <- function(x) {
    degpd_density(x, kappa = kappa, sigma = sigma, xi = xi)
  }

  # Joint density of (R, Z), evaluated at R = r, Z = z
  # Note: X = R w(z), Y = R (1-w(z)) = R w(-z) where w(z) = plogis(z)
  g <- function(z, r) {

    # Stable logistic transformations
    w <- plogis(z)
    v <- plogis(-z)

    x <- r * w
    y <- r * v

    u1 <- Fm(x)
    u2 <- Fm(y)

    # Avoid evaluating the copula density exactly at its boundaries
    eps <- 1e-20
    u1 <- pmin(pmax(u1, eps), 1 - eps)
    u2 <- pmin(pmax(u2, eps), 1 - eps)

    cop_density <- dCopula(
      cbind(u1, u2),
      copula = cop
    )

    # Jacobian of (R, Z) -> (X, Y)
    jacobian <- r * w * v

    cop_density * fm(x) * fm(y) * jacobian
  }

  # Symmetry: integrate over z >= 0 only.
  # The factor of 2 cancels between numerator and denominator.
  vapply(r, function(ri) {

    denominator <- integrate(
      f = function(z) g(z, ri),
      lower = -Inf,
      upper = Inf,
      rel.tol = rel.tol,
      subdivisions = subdivisions,
      stop.on.error = TRUE
    )$value

    numerator <- integrate(
      f = function(z) z^2 * g(z, ri),
      lower = -Inf,
      upper = Inf,
      rel.tol = rel.tol,
      subdivisions = subdivisions,
      stop.on.error = TRUE
    )$value

    if (!is.finite(denominator) || denominator <= 0) {
      stop("Invalid normalizing integral at radius ", ri)
    }

    # print(paste0("Numerator: ", numerator))
    # print(paste0("Denominator: ", denominator))

    sqrt(numerator / denominator)

  }, numeric(1))
}
