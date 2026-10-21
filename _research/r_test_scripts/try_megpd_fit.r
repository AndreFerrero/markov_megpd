library(copula)
source("megpd/code/megpd_fit.r")

gumbelC <- gumbelCopula(param = 5, dim = 2)
gaussC <- normalCopula(param = 0.8, dim = 2)
n <- 1e7

set.seed(46)
u <- rCopula(n, gaussC)
x <- egpd::qegpd(u[,1], kappa = 2, sigma = 1, xi = 0.1)
y <- egpd::qegpd(u[,2], kappa = 2, sigma = 1, xi = 0.1)

fit_results <- fit_megpd_bivariate(x, y, k = 10)

fit_radial <- fit_results$radial_fit
fit_delta <- fit_results$angular_fit

summary(fit_radial)
# plot(fit_radial)

summary(fit_delta)

# plot(fit_delta)

R <- x + y

V <- log(x / y)

# plot(R, V,
#      pch = 16, col = rgb(0,0,0,0.3),
#      main = "Angular spread vs radial component",
#      xlab = "R", ylab = "V")


plot_delta(fit_delta, R)
