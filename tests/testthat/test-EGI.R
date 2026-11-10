# End-to-end sequential design with EGI / EGIparallel.
#
# With fixed integration points, the discrete optimiser and no covariance
# re-estimation everything is deterministic: the first batch must be the one
# chosen by the corresponding max_* function on the initial model.

skip_on_cran()

T <- 80
m <- model_2d()
lower <- c(0, 0)
upper <- c(1, 1)
ip <- test_points(m, 30, seed = 21)
cand <- test_points(m, 15, seed = 50)
oc <- list(method = "discrete", optim.points = cand)
ic <- list(integration.points = ip)
kc <- list(CovReEstimate = FALSE)
iparam <- list(integration.points = ip, integration.weights = NULL)

run_egi <- function(...) {
  out <- NULL
  capture.output(out <- EGIparallel(T = T, model = m, fun = branin, lower = lower, upper = upper,
                                    optimcontrol = oc, integcontrol = ic, kmcontrol = kc, ...))
  out
}

expect_valid_run <- function(out, n_new, label) {
  expect_equal(nrow(out$par), n_new, label = paste(label, "number of new points"))
  expect_true(all(out$par >= 0 & out$par <= 1), label = paste(label, "inside the domain"))
  expect_equal(as.numeric(out$value), branin_vec(out$par), label = paste(label, "values = fun(par)"))
  expect_equal(out$lastmodel@n, m@n + n_new, label = paste(label, "model size"))
  expect_equal(unname(out$lastmodel@X[(m@n + 1):(m@n + n_new), , drop = FALSE]), unname(as.matrix(out$par)),
               label = paste(label, "new points in the model"))
}
same_par <- function(a, b) expect_equal(unname(as.matrix(a)), unname(as.matrix(b)), ignore_attr = TRUE)

## ---- batch-sequential criteria -------------------------------------------------------

test_that("EGIparallel with sur, jn, timse, imse and vorob starts with the max_* batch", {
  first <- list(
    sur = max_sur_parallel(lower, upper, oc, batchsize = 2, integration.param = iparam, T = T, model = m),
    jn = max_sur_parallel(lower, upper, oc, batchsize = 2, integration.param = iparam, T = T, model = m,
                          real.volume.variance = TRUE),
    timse = max_timse_parallel(lower, upper, oc, batchsize = 2, integration.param = iparam, T = T, model = m),
    imse = max_timse_parallel(lower, upper, oc, batchsize = 2, integration.param = iparam, T = T, model = m,
                              imse = TRUE),
    vorob = max_vorob_parallel(lower, upper, oc, batchsize = 2, integration.param = iparam, T = T, model = m))
  for (method in names(first)) {
    out <- run_egi(method = method, iter = 1, batchsize = 2)
    expect_valid_run(out, 2, method)
    same_par(out$par, first[[method]]$par)
    expect_equal(out$npoints, 2)
    expect_equal(out$nsteps, 2)
  }
})

test_that("EGIparallel runs several iterations", {
  out <- run_egi(method = "sur", iter = 2, batchsize = 2)
  expect_valid_run(out, 4, "sur, 2 iterations")
  expect_equal(out$nsteps, 4)
})

test_that("EGIparallel uses sur by default", {
  same_par(run_egi(iter = 1, batchsize = 2)$par, run_egi(method = "sur", iter = 1, batchsize = 2)$par)
})

test_that("EGIparallel with discrete candidate batches (batchsize * d columns)", {
  set.seed(51)
  batches <- matrix(runif(10 * 2 * 2), ncol = 4)
  first <- max_sur_parallel(lower, upper, list(method = "discrete", optim.points = batches), batchsize = 2,
                            integration.param = iparam, T = T, model = m)
  capture.output(out <- EGIparallel(T = T, model = m, method = "sur", fun = branin, iter = 1, batchsize = 2,
                                    lower = lower, upper = upper, integcontrol = ic, kmcontrol = kc,
                                    optimcontrol = list(method = "discrete", optim.points = batches)))
  expect_valid_run(out, 2, "full-batch candidates")
  same_par(out$par, first$par)

  first_v <- max_vorob_parallel(lower, upper, list(method = "discrete", optim.points = batches), batchsize = 2,
                                integration.param = iparam, T = T, model = m)
  capture.output(out_v <- EGIparallel(T = T, model = m, method = "vorob", fun = branin, iter = 1, batchsize = 2,
                                      lower = lower, upper = upper, integcontrol = ic, kmcontrol = kc,
                                      optimcontrol = list(method = "discrete", optim.points = batches)))
  expect_valid_run(out_v, 2, "vorob, full-batch candidates")
  same_par(out_v$par, first_v$par)
})

test_that("conservative criteria with a given Vorob'ev level start with the max_* batch", {
  mp <- list(penalization = 1, typeEx = ">", consLevel = 0.95, n_discrete_design = 50, consVorbLevel = 0.9)
  first_cons <- max_vorob_parallel(lower, upper, oc, batchsize = 2, integration.param = c(iparam, alpha = 0.9),
                                   T = T, model = m, penalisation = 1)
  first_vol <- max_futureVol_parallel(lower, upper, oc, batchsize = 2,
                                      integration.param = c(iparam, alpha = 0.9), T = T, model = m)
  set.seed(2)
  out_cons <- run_egi(method = "vorobCons", iter = 1, batchsize = 2, method.param = mp)
  set.seed(2)
  out_vol <- run_egi(method = "vorobVol", iter = 1, batchsize = 2, method.param = mp)
  expect_valid_run(out_cons, 2, "vorobCons")
  expect_valid_run(out_vol, 2, "vorobVol")
  same_par(out_cons$par, first_cons$par)
  same_par(out_vol$par, first_vol$par)
  expect_equal(out_cons$allCE_lvs[1], 0.9)
})

test_that("conservative criteria compute the conservative level at each iteration", {
  mp <- list(penalization = 1, typeEx = ">", consLevel = 0.95, n_discrete_design = 50)
  for (method in c("vorobCons", "vorobVol")) {
    set.seed(3)
    out <- run_egi(method = method, iter = 2, batchsize = 1, method.param = mp)
    expect_valid_run(out, 2, method)
    expect_length(out$allCE_lvs, 3)
    expect_true(all(is.finite(out$allCE_lvs) & out$allCE_lvs > 0 & out$allCE_lvs <= 1), label = method)
  }
})

test_that("criteria not available with several thresholds return 0 with a message", {
  for (method in c("vorob", "vorobCons", "vorobVol")) {
    expect_output(out <- EGIparallel(T = c(40, 80), model = m, method = method, fun = branin, iter = 1,
                                     batchsize = 1, lower = lower, upper = upper, optimcontrol = oc,
                                     integcontrol = ic, kmcontrol = kc), "not available with multiple thresholds")
    expect_equal(out, 0, label = method)
  }
})

test_that("sur and timse work with two thresholds", {
  for (method in c("sur", "timse")) {
    capture.output(out <- EGIparallel(T = c(40, 80), model = m, method = method, fun = branin, iter = 1,
                                      batchsize = 2, lower = lower, upper = upper, optimcontrol = oc,
                                      integcontrol = ic, kmcontrol = kc))
    expect_valid_run(out, 2, paste(method, "two thresholds"))
  }
})

test_that("noisy observations are added with new.noise.var", {
  mn <- model_2d(noise.var = 100)
  capture.output(out <- EGIparallel(T = T, model = mn, method = "sur", fun = branin, iter = 1, batchsize = 2,
                                    lower = lower, upper = upper, new.noise.var = 50, optimcontrol = oc,
                                    integcontrol = ic, kmcontrol = kc))
  expect_equal(out$lastmodel@n, mn@n + 2)
  expect_equal(tail(out$lastmodel@noise.var, 2), c(50, 50))
})

## ---- pointwise criteria (EGI) ------------------------------------------------------------

run_egi1 <- function(...) {
  out <- NULL
  capture.output(out <- EGI(T = T, model = m, fun = branin, lower = lower, upper = upper,
                            optimcontrol = oc, kmcontrol = kc, ...))
  out
}

test_that("EGI with ranjan, bichon and tmse starts with the max_infill_criterion point", {
  for (method in c("ranjan", "bichon", "tmse")) {
    first <- max_infill_criterion(lower, upper, oc, method = method, T = T, model = m)
    out <- run_egi1(method = method, iter = 2)
    expect_valid_run(out, 2, method)
    same_par(out$par[1, , drop = FALSE], first$par)
  }
})

test_that("EGI uses ranjan by default", {
  same_par(run_egi1(iter = 1)$par, run_egi1(method = "ranjan", iter = 1)$par)
})

test_that("EGI with tmse works with two thresholds, other criteria return 0", {
  capture.output(out <- EGI(T = c(40, 80), model = m, method = "tmse", fun = branin, iter = 1,
                            lower = lower, upper = upper, optimcontrol = oc, kmcontrol = kc))
  expect_valid_run(out, 1, "tmse two thresholds")
  expect_output(out0 <- EGI(T = c(40, 80), model = m, method = "ranjan", fun = branin, iter = 1,
                            lower = lower, upper = upper, optimcontrol = oc, kmcontrol = kc),
                "Multiple thresholds are supported only")
  expect_equal(out0, 0)
})

test_that("EGI with the tsee criterion", {
  first <- max_infill_criterion(lower, upper, oc, method = "tsee", T = T, model = m)
  out <- run_egi1(method = "tsee", iter = 1)
  expect_valid_run(out, 1, "tsee")
  same_par(out$par, first$par)
})

test_that("EGIparallel with a pointwise criterion switches to EGI with batchsize 1", {
  expect_output(out <- EGIparallel(T = T, model = m, method = "ranjan", fun = branin, iter = 1, batchsize = 2,
                                   lower = lower, upper = upper, optimcontrol = oc, kmcontrol = kc),
                "Switching to bachsize = 1")
  expect_equal(out$npoints, 1)
  same_par(out$par, run_egi1(method = "ranjan", iter = 1)$par)
})

test_that("EGI and EGIparallel dispatch to each other", {
  # EGI with a batch criterion calls EGIparallel with batchsize 1.
  capture.output(out <- EGI(T = T, model = m, method = "sur", fun = branin, iter = 1, lower = lower,
                            upper = upper, optimcontrol = oc, integcontrol = ic, kmcontrol = kc))
  same_par(out$par, run_egi(method = "sur", iter = 1, batchsize = 1)$par)
})
