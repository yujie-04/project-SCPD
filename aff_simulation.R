rm(list = ls())
set.seed(3)

case <- "drift"   # options: "jump", "two_jumps", "drift", "jump_var"

n <- 800
sigma <- 1

true_mean <- rep(0, n)
cp  <- NA
cp1 <- NA
cp2 <- NA

if (case == "jump") {
  cp <- 400
  true_mean <- c(rep(0, cp), rep(10, n - cp))
  x <- rnorm(n, mean = true_mean, sd = sigma)
  
} else if (case == "two_jumps") {
  cp1 <- 250
  cp2 <- 550
  true_mean <- c(rep(0, cp1), rep(10, cp2 - cp1), rep(3, n - cp2))
  x <- rnorm(n, mean = true_mean, sd = sigma)
  
} else if (case == "drift") {
  true_mean <- seq(0, 10, length.out = n)
  x <- rnorm(n, mean = true_mean, sd = sigma)
  
} else if (case == "jump_var") {
  cp <- 400
  true_mean <- c(rep(0, cp), rep(10, n - cp))
  sigma_t <- c(rep(1, cp), rep(2, n - cp))
  x <- rnorm(n, mean = true_mean, sd = sigma_t)
  
} else {
  stop("Unknown case. Use: jump, two_jumps, drift, jump_var")
}

eta <- 0.002        # learning rate
lambda0 <- 0.99    
lambda_min <- 0.60  

lambda <- rep(NA_real_, n)
xbar   <- rep(NA_real_, n)

m <- 0; w <- 0
Delta <- 0; Omega <- 0

# Initialise
lambda[1] <- lambda0
m <- x[1]
w <- 1
xbar[1] <- m / w

# AFF loop
for (t in 2:n) {
  xbar_prev <- xbar[t - 1]
  
  # keep previous states
  m_prev <- m; w_prev <- w
  Delta_prev <- Delta; Omega_prev <- Omega
  
  # L_t
  d_xbar <- (Delta_prev - xbar_prev * Omega_prev) / w_prev
  g <- 2 * (xbar_prev - x[t]) * d_xbar    
  
  # update lambda 
  lam_new <- lambda[t - 1] - eta * g
  
  # clamp to [lambda_min, 1]
  lam_new <- max(lambda_min, min(1, lam_new))
  lambda[t] <- lam_new
  
  # update mean states
  m <- lam_new * m_prev + x[t]
  w <- lam_new * w_prev + 1
  xbar[t] <- m / w
  
  # update derivative states
  Delta <- lam_new * Delta_prev + m_prev
  Omega <- lam_new * Omega_prev + w_prev
}


# Plot 
op <- par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

# mean tracking
plot(x, type = "l", col = "grey60", lwd = 1,
     main = "Adaptive Forgetting Factor Mean Tracking",
     xlab = "Time", ylab = "Value")
lines(true_mean, col = "blue", lwd = 2)
lines(xbar, col = "orange", lwd = 2, lty = 2)

# draw change points only when they exist
if (case == "jump" || case == "jump_var") abline(v = cp,  lty = 3, lwd = 2)
if (case == "two_jumps") {
  abline(v = cp1, lty = 3, lwd = 2)
  abline(v = cp2, lty = 3, lwd = 2)
}

# legend
legend_items <- c("Observed data", "True mean", "AFF estimate")
legend_cols  <- c("grey60", "blue", "orange")
legend_lty   <- c(1, 1, 2)
legend_lwd   <- c(1, 2, 2)

if (case != "drift") {
  legend_items <- c(legend_items, "Change points")
  legend_cols  <- c(legend_cols, "black")
  legend_lty   <- c(legend_lty, 3)
  legend_lwd   <- c(legend_lwd, 1.5)
}

legend("topleft",
       legend = legend_items,
       col    = legend_cols,
       lty    = legend_lty,
       lwd    = legend_lwd,
       bty    = "n",
       cex    = 0.65,
       seg.len = 1.5,
       x.intersp = 0.5,
       y.intersp = 0.5,
       inset  = c(0.01, 0.01))

# lambda
plot(lambda, type = "l", col = "darkgreen", lwd = 2,
     main = expression(paste("Adaptive forgetting factor ", lambda[t])),
     xlab = "Time", ylab = expression(lambda[t]),
     ylim = c(lambda_min, 1))
abline(h = lambda_min, lty = 2, col = "grey50")

if (case == "jump" || case == "jump_var") abline(v = cp,  lty = 3, lwd = 2)
if (case == "two_jumps") {
  abline(v = cp1, lty = 3, lwd = 2)
  abline(v = cp2, lty = 3, lwd = 2)
}

par(op)
