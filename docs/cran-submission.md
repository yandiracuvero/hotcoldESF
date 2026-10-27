# Submitting hotcoldESF to CRAN

Run these in R, from the repository root (`~/Github/Hotspot`), on your Mac.

## 0. Once

```r
install.packages(c("devtools", "rhub", "spelling", "urlchecker"))
```

## 1. Local checks

```r
devtools::document()            # regenerate man/ and NAMESPACE
devtools::spell_check()         # should report no errors (inst/WORDLIST lists names and terms)
urlchecker::url_check()         # links and DOIs in the documentation
devtools::check(remote = TRUE, manual = TRUE)   # full R CMD check --as-cran, with the PDF manual
```

## 2. Checks on other platforms

```r
devtools::check_win_devel()     # results by e-mail to yandira.cuvero@gmail.com (about 30 min)
devtools::check_win_release()
```

R-hub v2 runs on GitHub Actions, so it needs the package in a GitHub repository:

```r
rhub::rhub_setup()              # once: adds .github/workflows/rhub.yaml; commit and push it
rhub::rhub_doctor()             # checks the setup
rhub::rhub_check(platforms = c("linux", "macos-arm64", "windows"))
```

Fill in the results in `cran-comments.md`.

## 3. Before submitting

- Add `URL:` and `BugReports:` to DESCRIPTION when the GitHub repository is public.
- Check the DOI of Barde et al. (10.1093/ectj/utaf008) in `?hot_cold_map`.
- Bump nothing: the first version is 0.1.0.

## 4. Submit

```r
devtools::release()             # asks a few questions, builds and uploads
```

or upload the file from `devtools::build()` at https://cran.r-project.org/submit.html.
CRAN then sends a confirmation link to the maintainer address; click it to complete the submission.
After acceptance, tag the commit (`git tag v0.1.0`) and add `inst/CITATION` once the paper is published.
