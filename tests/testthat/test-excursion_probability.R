mn <- c(-2, -0.5, 0, 0.3, 1.7, 4)
sn <- c(0.5, 1, 2, 0.1, 1.5, 3)

test_that("one threshold gives P(Y > T) for Y ~ N(mn, sn^2)", {
  expect_equal(excursion_probability(mn, sn, T = 0.2), pnorm((mn - 0.2) / sn))
  expect_equal(excursion_probability(mn, sn, T = 0.2),
               pnorm(0.2, mean = mn, sd = sn, lower.tail = FALSE))
})

test_that("two thresholds give P(T1 < Y < T2)", {
  expected <- pnorm(1, mn, sn) - pnorm(-1, mn, sn)
  expect_equal(excursion_probability(mn, sn, T = c(-1, 1)), expected)
})

test_that("three thresholds give P(T1 < Y < T2) + P(Y > T3)", {
  expected <- (pnorm(0, mn, sn) - pnorm(-1, mn, sn)) + pnorm(2, mn, sn, lower.tail = FALSE)
  expect_equal(excursion_probability(mn, sn, T = c(-1, 0, 2)), expected)
})

test_that("the order of the thresholds does not matter", {
  expect_equal(excursion_probability(mn, sn, T = c(1, -1)),
               excursion_probability(mn, sn, T = c(-1, 1)))
})

test_that("probabilities are in [0, 1]", {
  for (T in list(0, c(-1, 1), c(-1, 0, 2))) {
    p <- excursion_probability(mn, sn, T = T)
    expect_true(all(p >= 0 & p <= 1))
  }
})
