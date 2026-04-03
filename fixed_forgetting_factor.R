rm(list = ls())
set.seed(123)

# 100 training data

n <- 100

# Stable mean and variance
# x_phase1 <- rnorm(n, mean = 0, sd = 1)

# Random walk
# x_phase1 <- cumsum(0.05 + rnorm(n, sd = 0.3))

# Abrupt mean shift
# x_phase1 <- c(rnorm(n/2, 0, 1), rnorm(n/2, 3, 1))

# Weak linear trend
 x_phase1 <- 0.02*(1:n) + rnorm(n, 0, 0.1)


# Forgetting factor mean recursion (finite)

ff_mean <- function(x, lambda) {
  n <- length(x)
  xbar <- numeric(n)
  xbar[1] <- x[1]
  
  for (t in 2:n) {
    a_t <- lambda * (1 - lambda^(t-1)) / (1 - lambda^t)
    b_t <- (1 - lambda) / (1 - lambda^t)
    xbar[t] <- a_t * xbar[t-1] + b_t * x[t]
  }
  xbar
}

# Objective function Q1(lambda) 

Q1_lambda <- function(lambda, x) {
  xbar <- ff_mean(x, lambda)
  sum( (x[2:length(x)] - xbar[1:(length(x)-1)])^2 )
}


# Estimate lambda* by using L-BFGS-B

estimate_lambda_LBFGSB <- function(x, init = 0.9, eps = 1e-6) {
  obj <- function(lam) Q1_lambda(lam, x)
  
  fit <- optim(
    par    = init,
    fn     = obj,
    method = "L-BFGS-B",
    lower  = eps,
    upper  = 1 - eps
  )
  
  list(lambda_star = fit$par,                  # Estimated optimal lambda
       Q1_min      = fit$value,                # Minimum value of Obj
       conv        = fit$convergence,          # 0 success, other fail
       message     = fit$message)
}


# Estimate lambda* by BFGS 

estimate_lambda_BFGS <- function(x, init_lambda = 0.9) {
  
  init_theta <- log(init_lambda / (1 - init_lambda))
  
  obj_theta <- function(theta) {
    lam <- 1 / (1 + exp(-theta))    
    Q1_lambda(lam, x)
  }
  
  fit <- optim(
    par    = init_theta,
    fn     = obj_theta,
    method = "BFGS"
  )
  
  lam_star <- 1 / (1 + exp(-fit$par))
  
  list(lambda_star = lam_star,
       Q1_min      = fit$value,
       conv        = fit$convergence,
       message     = fit$message)
}

# Run both estimators

fit_lbfgs <- estimate_lambda_LBFGSB(x_phase1, init = 0.8)
fit_bfgs  <- estimate_lambda_BFGS(x_phase1, init_lambda = 0.8)

cat(" Data-driven fixed lambda estimation \n")
cat(sprintf("L-BFGS-B: lambda* = %.4f,  Q1_min = %.4f\n",
            fit_lbfgs$lambda_star, fit_lbfgs$Q1_min))
cat(sprintf("BFGS:     lambda* = %.4f,  Q1_min = %.4f\n\n",
            fit_bfgs$lambda_star, fit_bfgs$Q1_min))


# Visualize Q1(lambda) curve

grid <- seq(0.01, 0.999, length.out = 200)
Qvals <- sapply(grid, Q1_lambda, x = x_phase1)

plot(grid, Qvals, type="l", lwd=2,
     xlab=expression(lambda), ylab=expression(Q[1](lambda)),
     main="One step ahead Q1(lambda)")
abline(v = fit_lbfgs$lambda_star, lty=2, lwd=2)
text(fit_lbfgs$lambda_star, min(Qvals),
     labels = paste0("lambda*=", round(fit_lbfgs$lambda_star,3)),
     pos=4)
