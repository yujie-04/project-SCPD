rm(list = ls())

set.seed(2025)
suppressPackageStartupMessages(library(mgcv))

if (!exists("x_phase1")) {
  # x_phase1 <- rnorm(100, mean = 0, sd = 1)                                   # stable
  # x_phase1 <- rnorm(100, mean = seq(0, 3, length.out = 100), sd = 1)       # mean drift
   x_phase1 <- rnorm(100, mean = 0, sd = seq(0.6, 1.6, length.out = 100))     # variance drift
  # x_phase1 <- rnorm(100, mean = seq(0, 3, length.out = 100),
  #                       sd   = seq(0.6, 1.6, length.out = 100))              # both drift
}
stopifnot(length(x_phase1) == 100)
t_phase1 <- 1:100

# Phase I
mean_model <- gam(x_phase1 ~ s(t_phase1))               # Generalized Additive Model
                                                        # x_t = f(t) + ε_t
sum_mean   <- summary(mean_model)
p_mean     <- as.numeric(sum_mean$s.table[1, "p-value"])

res1       <- residuals(mean_model)
log_res2   <- log(res1^2 + 1e-8)                      #log((r_t)^2) ≈ log((σ_t)^2)
var_model  <- gam(log_res2 ~ s(t_phase1))
sum_var    <- summary(var_model)
p_var      <- as.numeric(sum_var$s.table[1, "p-value"])

cat("\nPhase I results: \n")
cat(sprintf("Mean stability:     p = %.4g\n", p_mean))
cat(sprintf("Variance stability: p = %.4g\n", p_var))
cat("Rule: p > 0.05 ⇒ stable; p ≤ 0.05 ⇒ unstable.\n")

mean_stable <- (p_mean > 0.05)
var_stable  <- (p_var  > 0.05)
do_phase2   <- mean_stable && var_stable

if (!do_phase2) {
  cat("\nConclusion: PHASE I is unstable, Phase II will not run.\n")
  if (!mean_stable) cat(" - Mean unstable \n")
  if (!var_stable)  cat(" - Variance unstable \n")
  # Show Phase I diagnostics only
  op <- par(mfrow = c(1,2))
  plot(mean_model, shade = TRUE, main = "Phase I: Smooth Mean μ(t)",
       xlab = "t (1..100)", ylab = "μ(t)")
  plot(var_model,  shade = TRUE, main = "Phase I: Smooth log-Var log(σ²(t))",
       xlab = "t (1..100)", ylab = "log(σ²(t))")
  par(op)
  cat("\n(Stop here. Diagnose special causes, refit baseline later.)\n")
  invisible(NULL)
} else {
  cat("\nConclusion: Phase I is stable, proceeding to Phase II.\n")
}

# Lock Phase I targets
mu_hat_phase1    <- mean(x_phase1)
sigma_hat_phase1 <- sd(x_phase1)
cat(sprintf("Locked targets: μ̂₀ = %.4f,  σ̂ = %.4f\n", mu_hat_phase1, sigma_hat_phase1))

# Phase I diagnostics even when stable
plot(mean_model, shade = TRUE, main = "Phase I: Smooth Mean μ(t)",
     xlab = "t (1..100)", ylab = "μ(t)")
plot(var_model,  shade = TRUE, main = "Phase I: Smooth log-Var log(σ²(t))",
     xlab = "t (1..100)", ylab = "log(σ²(t))")

# Phase II (EWMA monitoring)
if (!exists("x_phase2")) {
  # Let change point at 120
  x_phase2 <- c(
    rnorm(120, mean = mu_hat_phase1,           sd = sigma_hat_phase1),
    rnorm(200 - 120, mean = mu_hat_phase1 + 1, sd = sigma_hat_phase1)
  )
}
stopifnot(length(x_phase2) == 200)

lambda <- 0.2
L      <- 3
n2     <- length(x_phase2)

z <- numeric(n2)
z[1] <- mu_hat_phase1
for (t in 2:n2) {
  z[t] <- lambda * x_phase2[t] + (1 - lambda) * z[t - 1]
}

ucl <- lcl <- numeric(n2)
for (t in 1:n2) {
  s_zt   <- sigma_hat_phase1 * sqrt((lambda / (2 - lambda)) * (1 - (1 - lambda)^(2 * t)))
  ucl[t] <- mu_hat_phase1 + L * s_zt
  lcl[t] <- mu_hat_phase1 - L * s_zt
}

plot(1:n2, z, type = "l", lwd = 2, col = "blue",
     ylim = range(c(lcl, ucl, z)),
     main = "Phase II EWMA Monitoring (Targets from Phase I)",
     xlab = "Sample Number (Phase II)", ylab = "EWMA z_t")
lines(1:n2, ucl, col = "red", lwd = 1.5, lty = 2)
lines(1:n2, lcl, col = "red", lwd = 1.5, lty = 2)
abline(h = mu_hat_phase1, col = "darkgreen", lwd = 1.5)
legend("topleft",
       legend = c("EWMA z_t", "UCL/LCL", "Target Mean (Phase I)"),
       col = c("blue", "red", "darkgreen"),
       lty = c(1, 2, 1), lwd = c(2, 1.5, 1.5))

signals <- which(z > ucl | z < lcl)
if (length(signals) > 0) {
  cat("\nEWMA ALERT: out of control indices (first 15): ",
      paste(head(signals, 15), collapse = ", "),
      if (length(signals) > 15) " ..." else "", "\n", sep = "")
} else {
  cat("\nEWMA: No out of control signals in Phase II.\n")
}
