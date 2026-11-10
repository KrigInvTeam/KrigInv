# computeRealVolumeConstant returns sum_ij w_i w_j P(Y(x_i) > T, Y(x_j) > T),
# i.e. the second moment of the (weighted) excursion volume, E[vol^2].
# Reference: bivariate probabilities computed with mvtnorm::pmvnorm.

T <- 80
m <- model_2d()

joint_exceed <- function(mu, Sigma, i, j) {
  if (i == j) return(pnorm(T, mu[i], sqrt(Sigma[i, i]), lower.tail = FALSE))
  ij <- c(i, j)
  as.numeric(mvtnorm::pmvnorm(lower = c(T, T), upper = c(Inf, Inf),
                              mean = mu[ij], sigma = Sigma[ij, ij]))
}

test_that("unweighted case equals the mean of the pairwise joint exceedance probabilities", {
  ip <- test_points(m, 15, seed = 7)
  p <- predict_nobias_km(m, as.data.frame(ip), type = "UK", cov.compute = TRUE)
  mu <- as.numeric(p$mean)
  M <- nrow(ip)
  ref <- mean(outer(seq_len(M), seq_len(M), Vectorize(function(i, j) joint_exceed(mu, p$cov, i, j))))
  expect_equal(computeRealVolumeConstant(m, ip, NULL, T), ref, tolerance = 1e-5)
})

test_that("weighted case (one weight per point) uses products of weights", {
  ip <- test_points(m, 10, seed = 7)
  w <- (1:10) / sum(1:10)
  p <- predict_nobias_km(m, as.data.frame(ip), type = "UK", cov.compute = TRUE)
  mu <- as.numeric(p$mean)
  P <- outer(1:10, 1:10, Vectorize(function(i, j) joint_exceed(mu, p$cov, i, j)))
  expect_equal(computeRealVolumeConstant(m, ip, w, T), sum(outer(w, w) * P), tolerance = 1e-5)
})

test_that("paired case (2n points, n weights) sums weighted joint probabilities of the pairs", {
  ip <- test_points(m, 16, seed = 8)
  w <- rep(1 / 8, 8)
  p <- predict_nobias_km(m, as.data.frame(ip), type = "UK", cov.compute = TRUE)
  mu <- as.numeric(p$mean)
  ref <- sum(w * vapply(1:8, function(i) joint_exceed(mu, p$cov, i, i + 8), numeric(1)))
  expect_equal(computeRealVolumeConstant(m, ip, w, T), ref, tolerance = 1e-5)
})

test_that("the constant is the second moment of the excursion volume (Monte Carlo)", {
  ip <- test_points(m, 15, seed = 7)
  p <- predict_nobias_km(m, as.data.frame(ip), type = "UK", cov.compute = TRUE)
  set.seed(11)
  Y <- mvtnorm::rmvnorm(1e5, mean = as.numeric(p$mean), sigma = p$cov + diag(1e-8, 15))
  vol2 <- mean(rowMeans(Y > T)^2)
  expect_equal(computeRealVolumeConstant(m, ip, NULL, T), vol2, tolerance = 0.01)
})
