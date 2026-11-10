# KrigInv 1.4.3

## Bug fixes

* `EGI()` and `EGIparallel()` no longer emit spurious "NAs introduced by coercion" warnings with penalized `km` models.

* `tmse_optim()` no longer returns `NaN` at points where the kriging variance is zero.

* `vorob_threshold()` and `integration_design(distrib = "vorob")` now work when the excursion set is empty or covers the whole domain. Before, they returned an empty value or stopped with an error.

* The `"tsee"` criterion can now be used in `max_infill_criterion()` and `EGI()`.

* Discrete optimisation over candidate batches (`optimcontrol$optim.points` with `batchsize*d` columns) now works in `max_sur_parallel()`, `max_timse_parallel()`, `max_vorob_parallel()`, `max_futureVol_parallel()` and `EGIparallel()`.

* `print_uncertainty_1d()` and `print_uncertainty_2d()`: the conservative estimate (`consQuantile = list(consLevel = ...)`) now works with one threshold and uses `consQuantile$typeEx`. With several thresholds it is skipped with a message instead of failing.

## Documentation

* Updated references: the published version of Azzimonti et al. (2021) <doi:10.1080/00401706.2019.1693427>, DOIs for all journal articles, and corrected author names and titles.

* Examples previously in `\dontrun{}` are now in `\donttest{}` and are checked.

* Fixed Rd markup that caused a "Lost braces" NOTE in CRAN checks.

* The package help page (`?KrigInv`) now lists the current maintainer.

## Testing and infrastructure

* New test suite based on `testthat` (3rd edition).

* GitHub Actions for `R CMD check` and test coverage (coveralls).
