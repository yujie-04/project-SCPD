rm(list = ls())
set.seed(101)

# Simulated data with a mean shift
n <- 100
mu0 <- 0
mu1 <- 1.5
sigma <- 1
change_point <- 70

x <- c(
  rnorm(change_point, mu0, sigma),
  rnorm(n - change_point, mu1, sigma)
)

# EWMA parameters
lambda <- 0.2
L <- 3

# EWMA statistic
z <- numeric(n)
z[1] <- mu0

for (t in 2:n) {
  z[t] <- lambda * x[t] + (1 - lambda) * z[t - 1]
}

# Control limits
ucl <- numeric(n)
lcl <- numeric(n)

for (t in 1:n) {
  ucl[t] <- mu0 + L * sigma *
    sqrt((lambda / (2 - lambda)) * (1 - (1 - lambda)^(2 * t)))
  lcl[t] <- mu0 - L * sigma *
    sqrt((lambda / (2 - lambda)) * (1 - (1 - lambda)^(2 * t)))
}

# Plotting
plot(1:n, z, type = "l", lwd = 2, col = "blue",
     ylim = range(c(lcl, ucl, z)),
     main = "EWMA Control Chart",
     xlab = "Time", ylab = "EWMA Statistic")

lines(1:n, ucl, col = "red", lwd = 1.5, lty = 2)
lines(1:n, lcl, col = "red", lwd = 1.5, lty = 2)
abline(h = mu0, col = "darkgreen", lwd = 1.5)
abline(v = change_point, col = "orange", lwd = 2, lty = 3)

legend("topleft",
       legend = c("EWMA statistic",
                  "Upper/Lower control limits",
                  "Target mean",
                  "True change point"),
       col = c("blue", "red", "darkgreen", "orange"),
       lty = c(1, 2, 1, 3),
       lwd = c(2, 1.5, 1.5, 2),
       bty = "n",
       cex = 0.75,
       y.intersp = 0.8,
       x.intersp = 0.5,
       seg.len = 1.8)
