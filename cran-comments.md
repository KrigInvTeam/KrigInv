## Test Environments

* Local macOS 27.2 install, R 4.6.2
* GitHub Actions: macos-latest (release), windows-latest (release), ubuntu-latest (devel, release, oldrel-1)
* win-builder: R-devel (Windows Server 2022 x64)

## R CMD check results

0 errors | 0 warnings | 0 notes on all platforms except win-builder (R-devel), which gives 1 NOTE:

```
Possibly misspelled words in DESCRIPTION:
  Azzimonti (37:253)
  Bect (37:288)
  Ginsbourger (37:197, 37:264)
  Picheny (37:185)
  Richet (37:297)
```

These are the surnames of the authors of the cited references and are spelled correctly.

## Update of existing CRAN package

* This update fixes the "Lost braces" NOTE in the Rd files currently appearing in the CRAN checks.
* It also fixes several bugs, updates the references and adds a test suite (see NEWS.md).
* Examples previously in `\dontrun{}` are now in `\donttest{}`.

## Reverse dependencies

* We checked the 2 reverse dependencies (pGPx 0.1.5, profExtrema 0.2.2) against this version with `R CMD check --run-donttest`: no problems.
