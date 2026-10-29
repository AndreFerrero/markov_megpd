library(copula)
library(egpd)

# Choose copula and radius
cop <- huslerReissCopula(6)
r <- 1
xi_true <- 0.2

# Unnormalised density of (R, Z) at R = r
g <- function(z, xi = 0.2) {
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

# Normalising constant
C <- integrate(g, -Inf, Inf, xi_true)$value

# Conditional density of Z given R = r
f_cond <- function(z) g(z, xi_true) / C

# Plot against z
z <- seq(-15, 15, length.out=1000)
plot(z, f_cond(z), type="l", lwd=2,
     xlab="z", ylab=expression(f[Z|R](z~"|"~r)),
     main=paste("Conditional density of Z given R =", r))