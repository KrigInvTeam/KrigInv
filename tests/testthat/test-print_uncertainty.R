# print_uncertainty_* plot an uncertainty measure on a grid and return its mean.
# Reference: the same measure computed directly from predict_nobias_km on the
# same grid.

measure <- function(pred, T, type) {
  mk <- as.numeric(pred$mean)
  sk <- as.numeric(pred$sd)
  pn <- excursion_probability(mk, sk, T)
  switch(type,
    pn = pn,
    sur = pn * (1 - pn),
    timse = sk^2 * dnorm(T, mk, sk),
    imse = sk^2,
    vorob = {
      alpha <- vorob_threshold(pn)
      pn * (pn <= alpha) + (1 - pn) * (pn > alpha)
    })
}
types <- c("pn", "sur", "timse", "imse", "vorob")

## ---- 1d ---------------------------------------------------------------------------

m1 <- model_1d()
T1 <- 0
grid1 <- matrix(seq(0, 1, length = 100), ncol = 1)
pred1 <- predict_nobias_km(m1, grid1, type = "UK")

test_that("print_uncertainty_1d returns the grid mean of each measure", {
  for (type in types) {
    val <- with_null_device(print_uncertainty_1d(m1, T1, type = type, resolution = 100))
    expect_equal(val, mean(measure(pred1, T1, type)), label = type)
  }
})

test_that("print_uncertainty_1d plotting options do not change the result", {
  ref <- mean(measure(pred1, T1, "pn"))
  val <- with_null_device(print_uncertainty_1d(m1, T1, type = "pn", resolution = 100, new.points = 2,
                                               vorobmean = TRUE, krigmeanplot = TRUE, Tplot = TRUE,
                                               consQuantile = 0.9))
  expect_equal(val, ref)
})

test_that("print_uncertainty_1d with an unknown type falls back to pn", {
  expect_output(val <- with_null_device(print_uncertainty_1d(m1, T1, type = "foo", resolution = 100)),
                "unknown 'type' value")
  expect_equal(val, mean(measure(pred1, T1, "pn")))
})

# Reference conservative level: anMC::conservativeEstimate on the same grid, same seed.
ce_level <- function(model, design, T, typeEx) {
  pred <- predict_nobias_km(model, design, type = "UK", cov.compute = TRUE)
  pred$cov <- pred$cov + 1e-7 * diag(nrow(pred$cov))
  set.seed(1)
  anMC::conservativeEstimate(alpha = 0.95, pred = pred, design = design, threshold = T, pn = NULL,
                             type = typeEx, verb = 0, lightReturn = TRUE, algo = "GANMC")$lvs
}

test_that("print_uncertainty_1d returns the conservative level for both excursion types", {
  grid50 <- matrix(seq(0, 1, length = 50), ncol = 1)
  for (typeEx in c(">", "<")) {
    set.seed(1)
    val <- with_null_device(print_uncertainty_1d(m1, T1, type = "pn", resolution = 50,
                                                 consQuantile = list(consLevel = 0.95, typeEx = typeEx)))
    expect_named(val, c("integrated pn", "CEquantile"))
    expect_equal(val[["CEquantile"]], ce_level(m1, grid50, T1, typeEx), tolerance = 1e-6, label = typeEx)
  }
})

test_that("print_uncertainty_1d skips the conservative estimate with several thresholds", {
  expect_output(val <- with_null_device(print_uncertainty_1d(m1, c(-3, 0), type = "pn", resolution = 50,
                                                             consQuantile = list(consLevel = 0.95))),
                "not available with multiple thresholds")
  expect_length(val, 1)
})

## ---- 2d ---------------------------------------------------------------------------

m2 <- model_2d()
T2 <- 80
s <- seq(0, 1, length = 30)
pred2 <- predict_nobias_km(m2, expand.grid(s, s), type = "UK")

test_that("print_uncertainty_2d returns the grid mean of each measure", {
  for (type in types) {
    val <- with_null_device(print_uncertainty_2d(m2, T2, type = type, resolution = 30))
    expect_equal(val, mean(measure(pred2, T2, type)), label = type)
  }
})

test_that("print_uncertainty_2d plotting options do not change the result", {
  val <- with_null_device(print_uncertainty_2d(m2, T2, type = "pn", resolution = 30, new.points = 2,
                                               krigmeanplot = TRUE, vorobmean = TRUE, consQuantile = 0.9))
  expect_equal(val, mean(measure(pred2, T2, "pn")))
})

test_that("print_uncertainty_2d with an unknown type falls back to pn", {
  expect_output(val <- with_null_device(print_uncertainty_2d(m2, T2, type = "foo", resolution = 30)),
                "unknown value for the type argument")
  expect_equal(val, mean(measure(pred2, T2, "pn")))
})

test_that("print_uncertainty_2d returns the conservative level for both excursion types", {
  s15 <- seq(0, 1, length = 15)
  for (typeEx in c(">", "<")) {
    set.seed(1)
    val <- with_null_device(print_uncertainty_2d(m2, T2, type = "pn", resolution = 15,
                                                 consQuantile = list(consLevel = 0.95, typeEx = typeEx)))
    expect_named(val, c("integrated pn", "CEquantile"))
    expect_equal(val[["CEquantile"]], ce_level(m2, expand.grid(s15, s15), T2, typeEx), tolerance = 1e-6,
                 label = typeEx)
  }
})

test_that("print_uncertainty_2d skips the conservative estimate with several thresholds", {
  expect_output(val <- with_null_device(print_uncertainty_2d(m2, c(40, 80), type = "pn", resolution = 15,
                                                             consQuantile = list(consLevel = 0.95))),
                "not available with multiple thresholds")
  expect_length(val, 1)
})

## ---- nd ---------------------------------------------------------------------------

m3 <- model_3d()
T3 <- 4
res <- 8
nint <- 20

# Prediction points of print_uncertainty_nd for the pair (d1, d2): a res x res grid
# on (d1, d2) times nint Sobol points on the remaining coordinate. Returns the
# matrix of the measure (one column per grid cell).
nd_measure <- function(type, d1, d2) {
  z <- randtoolbox::sobol(nint, dim = 1)
  g <- expand.grid(z = z, a = seq(0, 1, length = res), b = seq(0, 1, length = res))
  P <- matrix(0, nrow(g), 3)
  P[, d1] <- g$a
  P[, d2] <- g$b
  P[, setdiff(1:3, c(d1, d2))] <- g$z
  colnames(P) <- colnames(m3@X)
  matrix(measure(predict_nobias_km(m3, P, type = "UK"), T3, type), nrow = nint)
}

test_that("print_uncertainty_nd (one pair) returns the mean of the projected measure", {
  for (type in c("pn", "sur", "timse", "imse")) {
    val <- with_null_device(print_uncertainty_nd(m3, T3, type = type, resolution = res,
                                                 nintegpoints = nint, pairs = c(1, 2)))
    expect_equal(val, mean(colMeans(nd_measure(type, 1, 2))), label = type)
  }
})

test_that("print_uncertainty_nd options 'max' and 'min' project with max and min", {
  M <- nd_measure("pn", 1, 3)
  for (opt in c("max", "min")) {
    val <- with_null_device(print_uncertainty_nd(m3, T3, type = "pn", resolution = res,
                                                 nintegpoints = nint, pairs = c(1, 3), option = opt))
    expect_equal(val, mean(apply(M, 2, opt)), label = opt)
  }
})

test_that("print_uncertainty_nd without pairs returns the value of the last pair", {
  val <- with_null_device(print_uncertainty_nd(m3, T3, type = "pn", resolution = res, nintegpoints = nint))
  expect_equal(val, mean(colMeans(nd_measure("pn", 2, 3))))
})

test_that("print_uncertainty_nd switches 'vorob' to 'pn' with a message", {
  expect_output(val <- with_null_device(print_uncertainty_nd(m3, T3, type = "vorob", resolution = res,
                                                             nintegpoints = nint, pairs = c(1, 2))),
                "not available in n dimensions")
  expect_equal(val, mean(colMeans(nd_measure("pn", 1, 2))))
})

test_that("print_uncertainty_nd refuses models in dimension 1 or 2", {
  expect_output(val <- with_null_device(print_uncertainty_nd(m2, T2)), "use print_uncertainty_2d")
  expect_equal(val, 0)
})

## ---- dispatcher -------------------------------------------------------------------------

test_that("print_uncertainty calls the 1d, 2d or nd version according to the dimension", {
  expect_equal(with_null_device(print_uncertainty(m1, T1, type = "sur", resolution = 100)),
               with_null_device(print_uncertainty_1d(m1, T1, type = "sur", resolution = 100)))
  expect_equal(with_null_device(print_uncertainty(m2, T2, type = "sur", resolution = 30)),
               with_null_device(print_uncertainty_2d(m2, T2, type = "sur", resolution = 30)))
  expect_equal(with_null_device(print_uncertainty(m3, T3, type = "sur", resolution = res,
                                                  nintegpoints = nint, pairs = c(1, 2))),
               with_null_device(print_uncertainty_nd(m3, T3, type = "sur", resolution = res,
                                                     nintegpoints = nint, pairs = c(1, 2))))
})
