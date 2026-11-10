# The pointwise sampling criteria are closed-form expectations under the
# kriging predictive distribution Y(x) ~ N(mk, sk^2). Each one is checked
# against a Monte Carlo estimate of the expectation it is defined as.

T <- 80
m <- model_2d()
x <- test_points(m, 6, seed = 5)
pred <- predict_nobias_km(m, as.data.frame(x), type = "UK")
mk <- as.numeric(pred$mean)
sk <- as.numeric(pred$sd)

# Monte Carlo estimate of E[h(Y, mk, sk)] at each point (common random numbers).
mc_expectation <- function(h, n = 2e5, seed = 10) {
  set.seed(seed)
  z <- rnorm(n)
  vapply(seq_along(mk), function(i) mean(h(mk[i] + sk[i] * z, mk[i], sk[i])), numeric(1))
}

test_that("the test points are informative (non-trivial sd and criteria)", {
  expect_true(all(sk > 1))
  expect_gt(sum(ranjan_optim(x, m, T) > 1), 2)
})

test_that("ranjan_optim is E[(alpha^2 sk^2 - (Y - T)^2)_+]", {
  alpha <- 1.5
  mc <- mc_expectation(function(Y, mu, s) pmax(alpha^2 * s^2 - (Y - T)^2, 0))
  expect_equal(as.numeric(ranjan_optim(x, m, T, method.param = alpha)), mc, tolerance = 0.005)
})

test_that("bichon_optim is E[(alpha sk - |Y - T|)_+]", {
  alpha <- 1.5
  mc <- mc_expectation(function(Y, mu, s) pmax(alpha * s - abs(Y - T), 0))
  expect_equal(as.numeric(bichon_optim(x, m, T, method.param = alpha)), mc, tolerance = 0.005)
})

test_that("tsee_optim is E[(Y - T)_+] * E[(T - Y)_+]", {
  up <- mc_expectation(function(Y, mu, s) pmax(Y - T, 0))
  down <- mc_expectation(function(Y, mu, s) pmax(T - Y, 0))
  expect_equal(as.numeric(tsee_optim(x, m, T)), up * down, tolerance = 0.02)
})

test_that("tmse_optim is sk^2 times the N(mk, sk^2 + eps^2) density at T", {
  for (eps in c(0, 5)) {
    expected <- sk^2 * dnorm(T, mean = mk, sd = sqrt(sk^2 + eps^2))
    expect_equal(as.numeric(tmse_optim(x, m, T, method.param = eps)), expected)
  }
})

test_that("tmse_optim with several thresholds sums the densities", {
  Ts <- c(40, 80)
  expected <- sk^2 * (dnorm(Ts[1], mk, sk) + dnorm(Ts[2], mk, sk))
  expect_equal(as.numeric(tmse_optim(x, m, Ts)), expected)
})

test_that("criteria vanish at design points of a noiseless model", {
  xd <- m@X[1:3, , drop = FALSE]
  expect_equal(as.numeric(ranjan_optim(xd, m, T)), rep(0, 3), tolerance = 1e-8)
  expect_equal(as.numeric(bichon_optim(xd, m, T)), rep(0, 3), tolerance = 1e-8)
  expect_equal(as.numeric(tsee_optim(xd, m, T)), rep(0, 3), tolerance = 1e-8)
})

test_that("tmse_optim vanishes at design points of a noiseless model", {
  xd <- m@X[1:3, , drop = FALSE]
  expect_equal(as.numeric(tmse_optim(xd, m, T)), rep(0, 3), tolerance = 1e-8)
})

test_that("criteria are non-negative", {
  xg <- test_points(m, 200, seed = 6)
  for (crit in list(ranjan_optim, bichon_optim, tsee_optim, tmse_optim)) {
    expect_true(all(crit(xg, m, T) >= 0))
  }
})

test_that("a matrix of points gives the same values as one point at a time", {
  for (crit in list(ranjan_optim, bichon_optim, tsee_optim)) {
    one_by_one <- vapply(seq_len(nrow(x)), function(i) as.numeric(crit(x[i, ], m, T)), numeric(1))
    expect_equal(as.numeric(crit(x, m, T)), one_by_one)
  }
})
