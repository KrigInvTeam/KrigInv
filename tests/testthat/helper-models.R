# Shared fixtures for the KrigInv tests.
#
# All models are built with fixed covariance parameters (coef.cov, coef.var), so
# no likelihood optimisation is run: the trend is estimated by GLS in closed form
# and every result is deterministic and platform independent.

branin_vec <- function(X) apply(as.matrix(X), 1, branin)

# 2d branin model, noiseless (default) or with homogeneous observation noise.
model_2d <- function(n = 12, formula = ~1, noise.var = NULL, seed = 9) {
  set.seed(seed)
  X <- data.frame(matrix(runif(2 * n), ncol = 2))
  y <- branin_vec(X)
  if (!is.null(noise.var)) noise.var <- rep(noise.var, n)
  km(formula, design = X, response = y, covtype = "matern3_2",
     coef.cov = c(0.5, 0.8), coef.var = 2e4, noise.var = noise.var,
     control = list(trace = FALSE))
}

# Refit `model` on its design augmented with (newX, newy), same covariance
# parameters; the trend is re-estimated by GLS. `new.noise.var` is the noise
# variance of the new observations (only for models with noise.var).
refit_augmented <- function(model, newX, newy, new.noise.var = NULL) {
  X <- rbind(as.data.frame(model@X), setNames(as.data.frame(newX), colnames(model@X)))
  noise <- NULL
  if (model@noise.flag) noise <- c(model@noise.var, rep(new.noise.var, nrow(newX)))
  km(model@trend.formula, design = X, response = c(model@y, newy),
     covtype = model@covariance@name, coef.cov = model@covariance@range.val,
     coef.var = model@covariance@sd2, noise.var = noise, control = list(trace = FALSE))
}

# Points in [0,1]^d with the column names of `model`.
test_points <- function(model, n, seed = 1) {
  set.seed(seed)
  P <- matrix(runif(n * model@d), ncol = model@d)
  colnames(P) <- colnames(model@X)
  P
}

# 1d model: f(x) = 10 x sin(6 x) on [0, 1], 8 points, fixed covariance parameters.
f_1d <- function(x) 10 * x * sin(6 * x)
model_1d <- function(n = 8, seed = 3) {
  set.seed(seed)
  X <- data.frame(x = sort(runif(n)))
  km(~1, design = X, response = f_1d(X$x), covtype = "matern5_2",
     coef.cov = 0.25, coef.var = 25, control = list(trace = FALSE))
}

# 3d model: f(x) = sum((4 x - 2)^2) on [0, 1]^3, 15 points, fixed covariance parameters.
f_3d <- function(X) rowSums((4 * as.matrix(X) - 2)^2)
model_3d <- function(n = 15, seed = 4) {
  set.seed(seed)
  X <- data.frame(matrix(runif(3 * n), ncol = 3))
  km(~1, design = X, response = f_3d(X), covtype = "matern5_2",
     coef.cov = c(0.6, 0.6, 0.6), coef.var = 30, control = list(trace = FALSE))
}

# Run code with a null graphics device (no Rplots.pdf left behind).
with_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}
