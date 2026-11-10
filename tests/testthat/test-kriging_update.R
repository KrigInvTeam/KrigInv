# The fast update chain precomputeUpdateData -> computeQuickKrigcov ->
# predict_update_km_parallel must give the same kriging mean and sd as a model
# refitted on the augmented design with the same covariance parameters
# (Chevalier and Ginsbourger, corrected kriging update formulae).

update_chain <- function(model, intpoints, newX, newy, new.noise.var = 0) {
  pre <- precomputeUpdateData(model, intpoints)
  old <- predict_nobias_km(model, as.data.frame(intpoints), type = "UK")
  pnew <- predict_nobias_km(model, as.data.frame(newX), type = "UK", cov.compute = TRUE)
  kn <- computeQuickKrigcov(model, intpoints, newX, pre, pnew$F.newdata, pnew$c)
  Sigma.r <- pnew$cov + diag(new.noise.var, nrow(newX))
  predict_update_km_parallel(newXmean = pnew$mean, newXvar = pnew$sd^2, newXvalue = newy,
                             Sigma.r = Sigma.r, newdata.oldmean = old$mean,
                             newdata.oldsd = old$sd, kn = kn)
}

expect_update_matches_refit <- function(model, q, new.noise.var = 0) {
  intpoints <- test_points(model, 25, seed = 3)
  newX <- test_points(model, q, seed = 4)
  newy <- branin_vec(newX)
  upd <- update_chain(model, intpoints, newX, newy, new.noise.var)
  ref <- predict_nobias_km(refit_augmented(model, newX, newy, new.noise.var),
                           as.data.frame(intpoints), type = "UK")
  expect_equal(as.numeric(upd$mean), as.numeric(ref$mean), tolerance = 1e-8)
  expect_equal(as.numeric(upd$sd), as.numeric(ref$sd), tolerance = 1e-8)
}

test_that("update of one point matches refit (constant trend)", {
  expect_update_matches_refit(model_2d(), q = 1)
})

test_that("update of a batch matches refit (constant trend)", {
  expect_update_matches_refit(model_2d(), q = 3)
})

test_that("update of a batch matches refit (linear trend)", {
  expect_update_matches_refit(model_2d(formula = ~.), q = 3)
})

test_that("update with noisy new observations matches refit", {
  expect_update_matches_refit(model_2d(noise.var = 100), q = 2, new.noise.var = 100)
})

test_that("the updated variance never exceeds the current one", {
  m <- model_2d()
  intpoints <- test_points(m, 25, seed = 3)
  newX <- test_points(m, 3, seed = 4)
  old <- predict_nobias_km(m, as.data.frame(intpoints), type = "UK")
  upd <- update_chain(m, intpoints, newX, branin_vec(newX))
  expect_true(all(upd$sd <= old$sd + 1e-10))
})

test_that("a non positive definite Sigma.r returns an error flag", {
  out <- predict_update_km_parallel(newXmean = 0, newXvar = 1, newXvalue = 0,
                                    Sigma.r = matrix(-1), newdata.oldmean = 0,
                                    newdata.oldsd = 1, kn = matrix(1))
  expect_identical(out, list(error = TRUE))
})
