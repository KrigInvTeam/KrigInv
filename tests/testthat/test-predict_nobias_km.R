# predict_nobias_km is DiceKriging's predict.km without the variance bias
# correction; it also returns the quantities needed by the update formulas.

compare_with_predict_km <- function(model, newdata, type) {
  p1 <- predict_nobias_km(model, newdata, type = type, se.compute = TRUE, cov.compute = TRUE)
  p2 <- predict.km(model, newdata, type = type, se.compute = TRUE, cov.compute = TRUE,
                   bias.correct = FALSE, checkNames = FALSE)
  expect_equal(as.numeric(p1$mean), as.numeric(p2$mean))
  expect_equal(as.numeric(p1$sd), as.numeric(p2$sd))
  expect_equal(unname(p1$cov), unname(p2$cov))
}

test_that("matches predict.km(bias.correct = FALSE), noiseless model", {
  m <- model_2d()
  newdata <- as.data.frame(test_points(m, 6))
  compare_with_predict_km(m, newdata, "UK")
  compare_with_predict_km(m, newdata, "SK")
})

test_that("matches predict.km(bias.correct = FALSE), linear trend", {
  m <- model_2d(formula = ~.)
  compare_with_predict_km(m, as.data.frame(test_points(m, 6)), "UK")
})

test_that("matches predict.km(bias.correct = FALSE), noisy model", {
  m <- model_2d(noise.var = 100)
  compare_with_predict_km(m, as.data.frame(test_points(m, 6)), "UK")
})

test_that("se.compute = FALSE returns the mean only", {
  m <- model_2d()
  newdata <- as.data.frame(test_points(m, 4))
  p <- predict_nobias_km(m, newdata, type = "UK", se.compute = FALSE)
  expect_equal(as.numeric(p$mean), as.numeric(predict_nobias_km(m, newdata, type = "UK")$mean))
  expect_null(p$sd)
})

test_that("newdata with the wrong number of columns is an error", {
  m <- model_2d()
  expect_error(predict_nobias_km(m, data.frame(x = 0.5), type = "UK"), "same numbers of columns")
})
