library(copula)
library(egpd)
library(here)

set.seed(46)
theta <- 4
cop <- gumbelCopula(theta)
xi_true <- 0.2

U <- rCopula(10000, cop)
X <- apply(U, 2, qegpd, kappa = 2, sigma = 1, xi = xi_true)
R <- X[, 1] + X[, 2]

library(ReIns)

xi_hill <- Hill(R)
xi_hill$gamma[1:100]

plot(xi_hill$k[1:100], xi_hill$gamma[1:100], type = "l")

library(evd)
xi_gpd <- fpot(R, quantile(R, 0.95))
xi_gpd$estimate[2]

confint(xi_gpd)

exi(R)
chiplot(R)
mrlplot(R)
plot(xi_gpd)
tcplot(R, tlim = c(quantile(R, 0.85), quantile(R,0.999)))
