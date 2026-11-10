# Optimisers of the sampling criteria (max_* functions).
#
# Discrete optimisation is deterministic, so its result is compared with an
# independent run of the same greedy algorithm on the criterion functions:
# at step j the candidate minimising (maximising for max_futureVol_parallel)
# the criterion of the batch (candidate, points chosen so far) is added.
# Genoud runs are only checked for shape, bounds and internal consistency.

T <- 80
m <- model_2d()
lower <- c(0, 0)
upper <- c(1, 1)
ip <- test_points(m, 30, seed = 21)
iparam <- list(integration.points = ip, integration.weights = NULL)
cand <- test_points(m, 15, seed = 50)
q <- 3

old <- predict_nobias_km(m, as.data.frame(ip), type = "UK")
oldmean <- as.numeric(old$mean)
oldsd <- as.numeric(old$sd)
precalc <- precomputeUpdateData(m, ip)

# Greedy batch selection with crit(batch) evaluated on batch = rbind(candidate, chosen).
greedy <- function(crit, cand, q, pick = which.min) {
  chosen <- NULL
  steps <- list()
  for (j in seq_len(q)) {
    vals <- vapply(seq_len(nrow(cand)), function(i) crit(rbind(cand[i, ], chosen)), numeric(1))
    chosen <- rbind(chosen, cand[pick(vals), ])
    steps[[j]] <- vals
  }
  list(par = unname(chosen), last = steps[[q]], steps = steps)
}
xvec <- function(X) as.numeric(t(X))
discrete <- function(points = cand) list(method = "discrete", optim.points = points)
expect_par <- function(out, ref) expect_equal(unname(as.matrix(out$par)), ref$par, ignore_attr = TRUE)

## ---- max_sur_parallel -------------------------------------------------------------

test_that("max_sur_parallel (sur) picks the greedy batch", {
  pn <- excursion_probability(oldmean, oldsd, T)
  crit <- function(X) sur_optim_parallel(xvec(X), ip, NULL, oldmean, oldsd, precalc, m, T,
                                         new.noise.var = 0, batchsize = nrow(X),
                                         current.sur = mean(pn * (1 - pn)))
  ref <- greedy(crit, cand, q)
  out <- max_sur_parallel(lower, upper, discrete(), batchsize = q, integration.param = iparam, T = T, model = m)
  expect_par(out, ref)
  expect_equal(as.numeric(out$value), min(ref$last))
  expect_equal(out$allvalues, ref$last)
})

test_that("max_sur_parallel with two thresholds picks the greedy batch", {
  Ts <- c(40, 80)
  pn <- excursion_probability(oldmean, oldsd, Ts)
  ai <- matrix(oldmean, nrow = 2, ncol = length(oldmean), byrow = TRUE) - Ts
  crit <- function(X) sur_optim_parallel(xvec(X), ip, NULL, oldmean, oldsd, precalc, m, Ts,
                                         new.noise.var = 0, batchsize = nrow(X),
                                         current.sur = mean(pn * (1 - pn)), ai_precalc = ai)
  ref <- greedy(crit, cand, 2)
  out <- max_sur_parallel(lower, upper, discrete(), batchsize = 2, integration.param = iparam, T = Ts, model = m)
  expect_par(out, ref)
})

test_that("max_sur_parallel (jn) picks the greedy batch and adds the volume constant", {
  ai <- matrix(oldmean - T, nrow = 1)
  crit <- function(X) jn_optim_parallel(xvec(X), ip, NULL, oldmean, oldsd, precalc, m, T,
                                        new.noise.var = 0, batchsize = nrow(X), current.sur = 0,
                                        ai_precalc = ai)
  ref <- greedy(crit, cand, 2)
  out <- max_sur_parallel(lower, upper, discrete(), batchsize = 2, integration.param = iparam,
                          T = T, model = m, real.volume.variance = TRUE)
  expect_par(out, ref)
  expect_equal(as.numeric(out$value), min(ref$last))

  oc <- c(discrete(), real.volume.constant = TRUE)
  out_c <- max_sur_parallel(lower, upper, oc, batchsize = 2, integration.param = iparam,
                            T = T, model = m, real.volume.variance = TRUE)
  const <- computeRealVolumeConstant(m, ip, NULL, T)
  expect_equal(as.numeric(out_c$value), min(ref$last) + const)
  expect_equal(out_c$allvalues, ref$last + const)
})

test_that("max_sur_parallel (jn) with paired 'jn' integration points picks the greedy batch", {
  set.seed(5)
  jparam <- integration_design(list(distrib = "jn", n.points = 20), d = 2, lower = lower,
                               upper = upper, model = m, T = T)
  P <- jparam$integration.points
  w <- jparam$integration.weights
  pj <- predict_nobias_km(m, as.data.frame(P), type = "UK")
  pre_j <- precomputeUpdateData(m, P)
  ai <- matrix(as.numeric(pj$mean) - T, nrow = 1)
  crit <- function(X) jn_optim_parallel(xvec(X), P, w, as.numeric(pj$mean), as.numeric(pj$sd), pre_j,
                                        m, T, new.noise.var = 0, batchsize = nrow(X),
                                        current.sur = 0, ai_precalc = ai)
  ref <- greedy(crit, cand, 2)
  out <- max_sur_parallel(lower, upper, discrete(), batchsize = 2, integration.param = jparam,
                          T = T, model = m, real.volume.variance = TRUE)
  expect_par(out, ref)
})

## ---- max_timse_parallel -----------------------------------------------------------

test_that("max_timse_parallel picks the greedy batch (timse with epsilon, and imse)", {
  for (case in list(list(eps = 5, imse = FALSE), list(eps = 0, imse = TRUE))) {
    weight <- if (case$imse) 1 else dnorm(T, oldmean, sqrt(oldsd^2 + case$eps^2))
    cur <- mean(weight * oldsd^2)
    crit <- function(X) timse_optim_parallel(xvec(X), ip, NULL, oldmean, oldsd, precalc, m, T,
                                             new.noise.var = 0, weight = weight,
                                             batchsize = nrow(X), current.timse = cur)
    ref <- greedy(crit, cand, q)
    out <- max_timse_parallel(lower, upper, discrete(), batchsize = q, integration.param = iparam,
                              T = T, model = m, epsilon = case$eps, imse = case$imse)
    expect_par(out, ref)
    expect_equal(as.numeric(out$value), min(ref$last))
  }
})

## ---- max_vorob_parallel -----------------------------------------------------------

test_that("max_vorob_parallel picks the greedy batch (both excursion types, penalisation)", {
  for (case in list(list(typeEx = ">", pen = 1), list(typeEx = "<", pen = 2))) {
    pn <- if (case$typeEx == ">") pnorm((oldmean - T) / oldsd) else pnorm((T - oldmean) / oldsd)
    alpha <- vorob_threshold(pn)
    cur <- mean(pn * (pn <= alpha) + case$pen * (1 - pn) * (pn > alpha))
    crit <- function(X) vorob_optim_parallel(xvec(X), ip, NULL, oldmean, oldsd, precalc, m, T,
                                             new.noise.var = 0, batchsize = nrow(X), alpha = alpha,
                                             current.vorob = cur, penalisation = case$pen,
                                             typeEx = case$typeEx)
    ref <- greedy(crit, cand, q)
    out <- max_vorob_parallel(lower, upper, discrete(), batchsize = q, integration.param = iparam,
                              T = T, model = m, penalisation = case$pen, typeEx = case$typeEx)
    expect_par(out, ref)
    expect_equal(as.numeric(out$value), min(ref$last))
  }
})

## ---- max_futureVol_parallel -------------------------------------------------------

test_that("max_futureVol_parallel picks the greedy batch maximising the future volume", {
  alpha <- 0.9
  crit <- function(X) vorobVol_optim_parallel(xvec(X), ip, NULL, oldmean, oldsd, precalc, m, T,
                                              new.noise.var = 0, batchsize = nrow(X), alpha = alpha,
                                              current.crit = 0)
  ref <- greedy(crit, cand, q, pick = which.max)
  out <- max_futureVol_parallel(lower, upper, discrete(), batchsize = q,
                                integration.param = c(iparam, alpha = alpha), T = T, model = m)
  expect_par(out, ref)
  expect_equal(as.numeric(out$value), max(ref$last))
  expect_equal(unname(out$allvalues), do.call(cbind, ref$steps))
})

## ---- discrete optimisation over whole batches (batchsize * d columns) --------------

set.seed(51)
batches <- matrix(runif(10 * 2 * 2), ncol = 4)   # 10 candidate batches of 2 points

test_that("discrete optimisation over full batches evaluates every candidate batch", {
  pn <- excursion_probability(oldmean, oldsd, T)
  crit <- function(xv) sur_optim_parallel(xv, ip, NULL, oldmean, oldsd, precalc, m, T,
                                          new.noise.var = 0, batchsize = 2,
                                          current.sur = mean(pn * (1 - pn)))
  vals <- apply(batches, 1, crit)
  out <- max_sur_parallel(lower, upper, discrete(batches), batchsize = 2, integration.param = iparam,
                          T = T, model = m)
  expect_equal(out$allvalues, vals)
  expect_equal(as.numeric(out$value), min(vals))
  expect_equal(matrix(as.numeric(unlist(out$par)), nrow = 2),
               matrix(batches[which.min(vals), ], nrow = 2, byrow = TRUE))
})

test_that("discrete optimisation over full batches returns a numeric batch matrix", {
  ctrl <- discrete(batches)
  outs <- list(
    sur = max_sur_parallel(lower, upper, ctrl, batchsize = 2, integration.param = iparam, T = T, model = m),
    timse = max_timse_parallel(lower, upper, ctrl, batchsize = 2, integration.param = iparam, T = T, model = m),
    futureVol = max_futureVol_parallel(lower, upper, ctrl, batchsize = 2,
                                       integration.param = c(iparam, alpha = 0.9), T = T, model = m))
  for (nm in names(outs)) {
    o <- outs[[nm]]
    pick <- if (nm == "futureVol") which.max else which.min
    expect_true(is.numeric(o$par), label = nm)
    expect_equal(unname(o$par), matrix(batches[pick(o$allvalues), ], nrow = 2, byrow = TRUE), label = nm)
  }
})

test_that("max_vorob_parallel discrete optimisation over full batches", {
  pn <- pnorm((oldmean - T) / oldsd)
  alpha <- vorob_threshold(pn)
  cur <- mean(pn * (pn <= alpha) + (1 - pn) * (pn > alpha))
  vals <- apply(batches, 1, function(xv)
    vorob_optim_parallel(xv, ip, NULL, oldmean, oldsd, precalc, m, T, new.noise.var = 0, batchsize = 2,
                         alpha = alpha, current.vorob = cur, penalisation = 1))
  out <- max_vorob_parallel(lower, upper, discrete(batches), batchsize = 2, integration.param = iparam,
                            T = T, model = m)
  expect_equal(out$allvalues, vals)
  expect_equal(as.numeric(out$value), min(vals))
  expect_true(is.numeric(out$par))
  expect_equal(unname(out$par), matrix(batches[which.min(vals), ], nrow = 2, byrow = TRUE))
})

## ---- max_infill_criterion -----------------------------------------------------------

test_that("max_infill_criterion (discrete) returns the maximiser of each pointwise criterion", {
  crits <- list(tmse = function(x) tmse_optim(x, m, T),
                ranjan = function(x) ranjan_optim(x, m, T, method.param = 1),
                bichon = function(x) bichon_optim(x, m, T, method.param = 1))
  for (method in names(crits)) {
    vals <- as.numeric(crits[[method]](as.data.frame(cand)))
    out <- max_infill_criterion(lower, upper, discrete(), method = method, T = T, model = m)
    expect_equal(unname(as.numeric(out$par)), unname(cand[which.max(vals), ]), label = method)
    expect_equal(as.numeric(out$value), max(vals), label = method)
    expect_equal(as.numeric(out$allvalues), vals, label = method)
  }
})

test_that("max_infill_criterion works with the 'tsee' criterion", {
  vals <- as.numeric(tsee_optim(as.data.frame(cand), m, T))
  out <- max_infill_criterion(lower, upper, discrete(), method = "tsee", T = T, model = m)
  expect_equal(unname(as.numeric(out$par)), unname(cand[which.max(vals), ]))
  expect_equal(as.numeric(out$value), max(vals))
})

## ---- genoud ------------------------------------------------------------------------------

# The tiny budget makes genoud warn that it hit the generation limit: expected here.
genoud_control <- function(...) list(method = "genoud", pop.size = 10, max.generations = 2,
                                     wait.generations = 1, BFGSburnin = 0, print.level = 0, ...)

test_that("genoud runs return batches of the right shape inside the domain", {
  suppressWarnings(capture.output({
    outs <- list(
      sur1 = max_sur_parallel(lower, upper, genoud_control(optim.option = 1), batchsize = 2,
                              integration.param = iparam, T = T, model = m),
      sur2 = max_sur_parallel(lower, upper, genoud_control(optim.option = 2), batchsize = 2,
                              integration.param = iparam, T = T, model = m),
      timse = max_timse_parallel(lower, upper, genoud_control(), batchsize = 2,
                                 integration.param = iparam, T = T, model = m),
      vorob = max_vorob_parallel(lower, upper, genoud_control(), batchsize = 2,
                                 integration.param = iparam, T = T, model = m),
      futureVol = max_futureVol_parallel(lower, upper, genoud_control(), batchsize = 2,
                                         integration.param = c(iparam, alpha = 0.9), T = T, model = m),
      ranjan = max_infill_criterion(lower, upper, genoud_control(), method = "ranjan", T = T, model = m)
    )
  }))
  for (nm in names(outs)) {
    o <- outs[[nm]]
    nb <- if (nm == "ranjan") 1 else 2
    expect_equal(dim(o$par), c(nb, 2), label = nm)
    expect_true(all(o$par >= 0 & o$par <= 1), label = nm)
    expect_true(is.finite(as.numeric(o$value)), label = nm)
  }
})

test_that("the genoud value is the criterion at the returned batch", {
  pn <- excursion_probability(oldmean, oldsd, T)
  suppressWarnings(capture.output(out <- max_sur_parallel(lower, upper, genoud_control(optim.option = 1),
                                         batchsize = 2, integration.param = iparam, T = T, model = m)))
  val <- sur_optim_parallel(xvec(out$par), ip, NULL, oldmean, oldsd, precalc, m, T,
                            new.noise.var = 0, batchsize = 2, current.sur = mean(pn * (1 - pn)))
  expect_equal(as.numeric(out$value), val)
})
