library(copula)
library(egpd)

# Choose copula and radius
cop <- gumbelCopula(6)
z <- 2
xi_true <- -0.2

# Joint density of (R, Z) at R = r
g <- function(r, xi = 0.2) {
  w <- plogis(z)
  v <- plogis(-z)

  x <- r * w
  y <- r * v

  u1 <- pmin(pmax(pegpd(x, kappa=2, sigma=1, xi=xi), 1e-14), 1-1e-14)
  u2 <- pmin(pmax(pegpd(y, kappa=2, sigma=1, xi=xi), 1e-14), 1-1e-14)

  dCopula(cbind(u1, u2), copula=cop) *
    degpd_density(x, kappa=2, sigma=1, xi=xi) *
    degpd_density(y, kappa=2, sigma=1, xi=xi) *
    r * w * v
}

# Marginal density of Z
C <- integrate(g, 0, Inf, xi_true)$value

# Conditional density of R given Z
f_cond <- function(r) g(r, xi_true) / C

# Plot against z
r <- seq(0, 10, length.out=1000)
plot(r, f_cond(r), type="l", lwd=2,
     xlab="r", ylab=expression(f[R|Z](r~"|"~z)),
     main=paste("Conditional density of R given Z =", z))