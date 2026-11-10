test_that("NULL input returns 0.5", {
  expect_equal(vorob_threshold(NULL), 0.5)
})

test_that("a constant coverage function returns that constant", {
  expect_equal(vorob_threshold(rep(0.3, 50)), 0.3)
})

test_that("the threshold lies within the range of pn", {
  set.seed(1)
  pn <- runif(200)
  alpha <- vorob_threshold(pn)
  expect_gte(alpha, min(pn))
  expect_lte(alpha, max(pn))
})

test_that("the Vorob'ev quantile set {pn >= alpha} has volume close to mean(pn)", {
  # Defining property of the Vorob'ev threshold, up to the discretisation 1/M.
  set.seed(2)
  for (shape in c(0.3, 1, 3)) {
    pn <- rbeta(500, shape, 1)
    alpha <- vorob_threshold(pn)
    expect_lte(abs(mean(pn >= alpha) - mean(pn)), 1 / length(pn) + 1e-12)
  }
})

test_that("a coverage function close to 1 everywhere returns a number", {
  expect_equal(vorob_threshold(rep(1, 10)), 1)
  expect_length(vorob_threshold(c(rep(1, 9), 0.99)), 1)
  pn <- c(rep(1, 9), 0.99)
  expect_equal(vorob_threshold(pn), min(pn))
  expect_lte(abs(mean(pn >= vorob_threshold(pn)) - mean(pn)), 1 / length(pn))
})
