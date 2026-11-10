# Regression test: EGI / EGIparallel used length(model@penalty == 0), which
# coerced the penalty list to numeric ("NAs introduced by coercion" warnings)
# for penalized models.

penalized_model <- function() {
  set.seed(1)
  X <- data.frame(matrix(runif(2 * 15), ncol = 2))
  km(~1, design = X, response = branin_vec(X), covtype = "gauss",
     penalty = list(fun = "SCAD", value = 3), control = list(trace = FALSE))
}

candidates <- as.matrix(expand.grid(seq(0.05, 0.95, length.out = 5),
                                    seq(0.05, 0.95, length.out = 5)))

test_that("EGI with a penalized model gives no warnings and keeps the penalty", {
  m <- penalized_model()
  expect_no_warning(capture.output(
    out <- EGI(T = 80, model = m, method = "ranjan", fun = branin, iter = 1,
               lower = c(0, 0), upper = c(1, 1),
               optimcontrol = list(method = "discrete", optim.points = candidates))))
  expect_equal(out$lastmodel@penalty$fun, "SCAD")
  expect_equal(out$lastmodel@penalty$value, 3)
})

test_that("EGIparallel with a penalized model gives no warnings and keeps the penalty", {
  m <- penalized_model()
  expect_no_warning(capture.output(
    out <- EGIparallel(T = 80, model = m, method = "sur", fun = branin, iter = 1,
                       batchsize = 2, lower = c(0, 0), upper = c(1, 1),
                       optimcontrol = list(method = "discrete", optim.points = candidates),
                       integcontrol = list(distrib = "sobol", n.points = 30))))
  expect_equal(out$lastmodel@penalty$fun, "SCAD")
  expect_equal(out$lastmodel@penalty$value, 3)
})
