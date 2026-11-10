T <- 80
m <- model_2d()
lower <- c(0, 0)
upper <- c(1, 1)
in_bounds <- function(P, lower, upper) all(t(P) >= lower & t(P) <= upper)

# Candidate set used as the base of the importance sampling ("spec" init.distrib).
cand <- test_points(m, 400, seed = 40)
pred_cand <- predict_nobias_km(m, as.data.frame(cand), type = "UK")
pn_cand <- excursion_probability(pred_cand$mean, pred_cand$sd, T)
sd_at <- function(P) as.numeric(predict_nobias_km(m, as.data.frame(P), type = "UK")$sd)

is_design <- function(distrib, n.points = 4000, min.prob = 0.5, seed = 1) {
  set.seed(seed)
  integration_design(list(distrib = distrib, n.points = n.points, n.candidates = nrow(cand),
                          init.distrib = "spec", init.distrib.spec = cand),
                     d = 2, lower = lower, upper = upper, model = m, T = T, min.prob = min.prob)
}

## ---- simple designs -----------------------------------------------------------

test_that("NULL integcontrol gives d * 100 Sobol points in the domain", {
  r <- integration_design(NULL, d = 2, lower = c(-1, 2), upper = c(3, 5), model = m)
  expect_equal(dim(r$integration.points), c(200, 2))
  expect_true(in_bounds(r$integration.points, c(-1, 2), c(3, 5)))
  expect_null(r$integration.weights)
  expect_equal(colnames(r$integration.points), colnames(m@X))
})

test_that("'sobol' gives the scaled Sobol sequence", {
  r <- integration_design(list(distrib = "sobol", n.points = 50), d = 2, lower = c(-1, 2), upper = c(3, 5))
  expected <- t(c(-1, 2) + t(randtoolbox::sobol(50, dim = 2)) * c(4, 3))
  expect_equal(unname(r$integration.points), unname(expected))
})

test_that("'MC' gives n.points uniform points in the domain", {
  set.seed(3)
  r <- integration_design(list(distrib = "MC", n.points = 70), d = 2, lower = c(-1, 2), upper = c(3, 5))
  expect_equal(dim(r$integration.points), c(70, 2))
  expect_true(in_bounds(r$integration.points, c(-1, 2), c(3, 5)))
})

test_that("pre-specified integration points and weights are returned unchanged", {
  P <- test_points(m, 10, seed = 41)
  w <- rep(0.1, 10)
  r <- integration_design(list(integration.points = P, integration.weights = w),
                          lower = lower, upper = upper, model = m)
  expect_equal(r$integration.points, P)
  expect_equal(r$integration.weights, w)
})

test_that("d = 1 returns a one-column matrix", {
  for (distrib in c("sobol", "MC")) {
    r <- integration_design(list(distrib = distrib, n.points = 20), d = 1, lower = 0, upper = 1)
    expect_true(is.matrix(r$integration.points))
    expect_equal(dim(r$integration.points), c(20, 1))
  }
})

## ---- input problems (reported with print() and NULL, not with stop()) --------

test_that("lower and upper of different lengths give NULL and a message", {
  expect_output(r <- integration_design(NULL, lower = c(0, 0), upper = 1), "must have the same length")
  expect_null(r)
})

test_that("importance sampling without a model gives NULL and a message", {
  for (distrib in c("sur", "vorob", "jn", "timse", "imse")) {
    expect_output(r <- integration_design(list(distrib = distrib, n.points = 10), d = 2,
                                          lower = lower, upper = upper, T = T),
                  "you must set the argument 'model'")
    expect_null(r)
  }
})

## ---- importance sampling designs ------------------------------------------------

test_that("importance sampling designs have the right size and positive weights", {
  for (distrib in c("sur", "vorob", "timse", "imse")) {
    r <- is_design(distrib, n.points = 50)
    expect_equal(dim(r$integration.points), c(50, 2))
    expect_length(r$integration.weights, 50)
    expect_true(all(r$integration.weights > 0))
    expect_equal(colnames(r$integration.points), colnames(m@X))
  }
})

test_that("'jn' returns 2n points (pairs) and n weights", {
  r <- is_design("jn", n.points = 50)
  expect_equal(dim(r$integration.points), c(100, 2))
  expect_length(r$integration.weights, 50)
})

test_that("with init.distrib = 'spec' all points are taken from the given candidates", {
  for (distrib in c("sur", "vorob", "jn", "timse", "imse")) {
    r <- is_design(distrib, n.points = 50)
    keys <- function(P) apply(signif(P, 12), 1, paste, collapse = ",")
    expect_true(all(keys(r$integration.points) %in% keys(cand)), label = distrib)
  }
})

test_that("init.distrib = 'MC' and 'sobol' also work", {
  for (init in c("MC", "sobol")) {
    set.seed(4)
    r <- integration_design(list(distrib = "sur", n.points = 30, init.distrib = init),
                            d = 2, lower = lower, upper = upper, model = m, T = T)
    expect_equal(dim(r$integration.points), c(30, 2))
    expect_true(in_bounds(r$integration.points, lower, upper))
  }
})

test_that("weighted sums are unbiased estimates of the average over the candidates", {
  # sum_i w_i f(x_i) estimates mean(f(candidates)). min.prob = 0.5 bounds the
  # weights, so the Monte Carlo error is small (about 1%).
  for (distrib in c("sur", "vorob", "timse", "imse")) {
    r <- is_design(distrib)
    expect_equal(sum(r$integration.weights), 1, tolerance = 0.04, label = paste(distrib, "sum of weights"))
    expect_equal(sum(r$integration.weights * sd_at(r$integration.points)), mean(pred_cand$sd),
                 tolerance = 0.04, label = paste(distrib, "weighted mean of sd"))
  }
})

test_that("'jn' pair weights estimate the product of two averages", {
  r <- is_design("jn")
  n <- length(r$integration.weights)
  P1 <- r$integration.points[1:n, ]
  P2 <- r$integration.points[(n + 1):(2 * n), ]
  expect_equal(sum(r$integration.weights * sd_at(P1) * sd_at(P2)), mean(pred_cand$sd)^2, tolerance = 0.04)
})

test_that("'sur' samples the candidates proportionally to pn (1 - pn)", {
  # Sampling probabilities prob_j proportional to max(tau_j / sum(tau), min.prob / N),
  # so the mean of tau over the sampled points estimates sum(prob_j * tau_j).
  min.prob <- 0.001
  tau <- pn_cand * (1 - pn_cand)
  prob <- pmax(tau / sum(tau), min.prob / length(tau))
  prob <- prob / sum(prob)
  r <- is_design("sur", n.points = 4000, min.prob = min.prob)
  p <- predict_nobias_km(m, as.data.frame(r$integration.points), type = "UK")
  pn_pts <- excursion_probability(p$mean, p$sd, T)
  expect_equal(mean(pn_pts * (1 - pn_pts)), sum(prob * tau), tolerance = 0.02)
  expect_gt(sum(prob * tau), mean(tau))
})

test_that("'vorob' returns the Vorob'ev threshold of the candidates", {
  r <- is_design("vorob", n.points = 50)
  expect_equal(r$alpha, vorob_threshold(pn_cand))
})

test_that("'vorob' works when the excursion set is (almost) empty or the whole domain", {
  for (Tx in c(-1000, 1000)) {
    set.seed(1)
    r <- integration_design(list(distrib = "vorob", n.points = 10, init.distrib = "spec",
                                 init.distrib.spec = cand[1:50, ]),
                            d = 2, lower = lower, upper = upper, model = m, T = Tx)
    p <- predict_nobias_km(m, as.data.frame(cand[1:50, ]), type = "UK")
    expect_equal(r$alpha, vorob_threshold(excursion_probability(p$mean, p$sd, Tx)))
  }
})
