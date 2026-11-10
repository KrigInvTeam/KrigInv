# Batch (parallel) SUR-type criteria.
#
# Reference: the batch observations Y ~ N(mk(X), Sigma(X)) are simulated jointly,
# the posterior mean is updated with the exact linear kriging update (checked
# against refits in test-kriging_update.R) and the future excursion probability
# p+ is computed. Each criterion must equal the Monte Carlo average of the
# quantity it is defined as:
#   sur      : E[ mean p+ (1 - p+) ]
#   jn       : - E[ (mean p+)^2 ]        (the constant E[vol^2] is added elsewhere)
#   vorob    : E[ mean( p+ 1{p+ < alpha} + pen (1 - p+) 1{p+ >= alpha} ) ]
#   vorobVol : E[ mean 1{p+ >= alpha} ]
#   timse    : mean( weight * sd+^2 )     (deterministic, checked against a refit)

T <- 80
m <- model_2d()
ip <- test_points(m, 30, seed = 21)
X <- test_points(m, 2, seed = 22)          # a batch of two new points
xvec <- function(X) as.numeric(t(X))       # criteria take the batch as one vector

old <- predict_nobias_km(m, as.data.frame(ip), type = "UK")
oldmean <- as.numeric(old$mean)
oldsd <- as.numeric(old$sd)
precalc <- precomputeUpdateData(m, ip)
pn <- excursion_probability(oldmean, oldsd, T)
alpha <- vorob_threshold(pn)

# Simulated future excursion probabilities: an M x N matrix (one column per draw).
future_p <- function(X, T, N = 4e4, seed = 30, typeEx = ">") {
  pb <- predict_nobias_km(m, as.data.frame(X), type = "UK", cov.compute = TRUE)
  kn <- computeQuickKrigcov(m, ip, X, precalc, pb$F.newdata, pb$c)
  upd <- predict_update_km_parallel(pb$mean, pb$sd^2, pb$mean, pb$cov, oldmean, oldsd, kn)
  set.seed(seed)
  Y <- mvtnorm::rmvnorm(N, mean = as.numeric(pb$mean), sigma = pb$cov, method = "chol")
  newmean <- oldmean + upd$lambda %*% (t(Y) - as.numeric(pb$mean))
  p <- excursion_probability(newmean, as.numeric(upd$sd), T)
  if (typeEx == "<") p <- 1 - p
  p
}

sur <- function(X, T = 80, weights = NULL, new.noise.var = NULL, model = m) {
  ai <- if (length(T) > 1) matrix(oldmean, nrow = length(T), ncol = length(oldmean), byrow = TRUE) - T else NULL
  sur_optim_parallel(xvec(X), ip, weights, oldmean, oldsd, precalc, model, T,
                     new.noise.var = new.noise.var, batchsize = nrow(X),
                     current.sur = mean(pn * (1 - pn)), ai_precalc = ai)
}
jn <- function(X, weights = NULL) {
  jn_optim_parallel(xvec(X), ip, weights, oldmean, oldsd, precalc, m, T, batchsize = nrow(X),
                    current.sur = 0, ai_precalc = matrix(oldmean - T, nrow = 1))
}
vorob <- function(X, typeEx = ">", pen = 1, weights = NULL) {
  vorob_optim_parallel(xvec(X), ip, weights, oldmean, oldsd, precalc, m, T, batchsize = nrow(X),
                       alpha = alpha, current.vorob = 1, penalisation = pen, typeEx = typeEx)
}
vorobVol <- function(X, typeEx = ">", weights = NULL) {
  vorobVol_optim_parallel(xvec(X), ip, weights, oldmean, oldsd, precalc, m, T, batchsize = nrow(X),
                          alpha = alpha, current.crit = 0, typeEx = typeEx)
}
timse_weight <- dnorm(T, oldmean, oldsd)
current_timse <- mean(timse_weight * oldsd^2)
timse <- function(X, weights = NULL, new.noise.var = 0) {
  timse_optim_parallel(xvec(X), ip, weights, oldmean, oldsd, precalc, m, T,
                       new.noise.var = new.noise.var, weight = timse_weight,
                       batchsize = nrow(X), current.timse = current_timse)
}

## ---- agreement with the Monte Carlo / refit references ----------------------

test_that("sur is E[mean p+(1 - p+)]", {
  p <- future_p(X, T)
  expect_equal(sur(X), mean(p * (1 - p)), tolerance = 0.003)
})

test_that("sur with two thresholds is E[mean p+(1 - p+)] for P(T1 < Y < T2)", {
  Ts <- c(40, 80)
  p <- future_p(X, Ts)
  expect_equal(sur(X, T = Ts), mean(p * (1 - p)), tolerance = 0.003)
})

test_that("jn is -E[(mean p+)^2]", {
  p <- future_p(X, T)
  expect_equal(jn(X), -mean(colMeans(p)^2), tolerance = 0.003)
})

test_that("vorob is the expected type II + penalised type I error", {
  for (typeEx in c(">", "<")) {
    for (pen in c(0, 1, 2)) {
      p <- future_p(X, T, typeEx = typeEx)
      ref <- mean(p * (p < alpha) + pen * (1 - p) * (p >= alpha))
      expect_equal(vorob(X, typeEx = typeEx, pen = pen), ref, tolerance = 0.003,
                   label = paste("vorob, typeEx", typeEx, "pen", pen))
    }
  }
})

test_that("vorobVol is the expected volume of the future set {p+ >= alpha}", {
  for (typeEx in c(">", "<")) {
    p <- future_p(X, T, typeEx = typeEx)
    expect_equal(vorobVol(X, typeEx = typeEx), mean(p >= alpha), tolerance = 0.003,
                 label = paste("vorobVol, typeEx", typeEx))
  }
})

test_that("timse is the weighted mean of the updated variance", {
  newsd <- as.numeric(predict_nobias_km(refit_augmented(m, X, branin_vec(X)),
                                        as.data.frame(ip), type = "UK")$sd)
  expect_equal(timse(X), mean(timse_weight * newsd^2), tolerance = 1e-8)
})

## ---- properties that hold for any batch -------------------------------------

test_that("sur and timse never exceed the current uncertainty", {
  for (seed in 1:5) {
    Xs <- test_points(m, 2, seed = 100 + seed)
    expect_lte(sur(Xs), mean(pn * (1 - pn)))
    expect_lte(timse(Xs), current_timse)
  }
})

test_that("adding a point to the batch does not increase sur, timse or jn", {
  X3 <- rbind(X, test_points(m, 1, seed = 23))
  expect_lte(sur(X3), sur(X))
  expect_lte(timse(X3), timse(X))
  expect_lte(jn(X3), jn(X))
})

test_that("uniform integration weights 1/M give the unweighted criteria", {
  w <- rep(1 / nrow(ip), nrow(ip))
  expect_equal(sur(X, weights = w), sur(X))
  expect_equal(timse(X, weights = w), timse(X))
  expect_equal(jn(X, weights = w), jn(X))
  expect_equal(vorob(X, weights = w), vorob(X))
  expect_equal(vorobVol(X, weights = w), vorobVol(X))
})

## ---- *_optim_parallel2 wrappers ---------------------------------------------

test_that("the *_optim_parallel2 wrappers equal the batch criteria", {
  x1 <- X[1, ]
  others <- X[2, ]
  ai1 <- matrix(oldmean - T, nrow = 1)
  expect_equal(
    sur_optim_parallel2(x1, others, ip, NULL, oldmean, oldsd, precalc, m, T, batchsize = 2,
                        current.sur = mean(pn * (1 - pn))), sur(X))
  expect_equal(
    jn_optim_parallel2(x1, others, ip, NULL, oldmean, oldsd, precalc, m, T, batchsize = 2,
                       current.sur = 0, ai_precalc = ai1), jn(X))
  expect_equal(
    timse_optim_parallel2(x1, others, ip, NULL, oldmean, oldsd, precalc, m, T,
                          new.noise.var = 0, weight = timse_weight, batchsize = 2,
                          current.timse = current_timse), timse(X))
  expect_equal(
    vorob_optim_parallel2(x1, others, ip, NULL, oldmean, oldsd, precalc, m, T, batchsize = 2,
                          alpha = alpha, current.vorob = 1, penalisation = 1), vorob(X))
  expect_equal(
    vorobVol_optim_parallel2(x1, others, ip, NULL, oldmean, oldsd, precalc, m, T, batchsize = 2,
                             alpha = alpha, current.crit = 0), vorobVol(X))
})

## ---- batch points on existing design points ----------------------------------

test_that("a batch point on a design point of a noiseless model returns the penalised current value", {
  Xd <- rbind(m@X[1, ], X[1, ])
  cur_sur <- mean(pn * (1 - pn))
  expect_equal(sur(Xd), cur_sur * 1.01)
  expect_equal(timse(Xd), current_timse * 1.01)
  expect_equal(jn(Xd), 0)              # current.sur passed as 0
  expect_equal(vorob(Xd), 1 + 0.01)    # current.vorob passed as 1
  expect_equal(vorobVol(Xd), 0 + 0.01) # current.crit passed as 0
})

## ---- noisy observations ------------------------------------------------------

test_that("new.noise.var does not enter the update used by the criteria (intended)", {
  # Design choice: the criteria update the model as if the new observations
  # were noiseless; new.noise.var only disables the duplicate-point check.
  expect_equal(sur(X, new.noise.var = 100), sur(X))
  expect_equal(timse(X, new.noise.var = 100), timse(X))
})

test_that("with new.noise.var, a design point is evaluated instead of penalised", {
  mn <- model_2d(noise.var = 100)
  Xd <- rbind(mn@X[1, ], X[1, ])
  val <- sur(Xd, new.noise.var = 100, model = mn)
  expect_true(is.finite(val))
  expect_false(isTRUE(all.equal(val, mean(pn * (1 - pn)) * 1.01)))
})
