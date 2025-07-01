## Test environments

* Local PC:
  - Windows 10, R devel

* [GitHub Actions](https://github.com/ms609/Rogue/actions)
  - ubuntu-latest
    - R 4.1
    - R release (tests, examples & vignettes run with valgrind)
    - R devel
  - macOS-latest, R release
  - windows-latest, R release
  - [R-hub](https://github.com/ms609/Rogue/actions/workflows/rhub.yaml)

* `devtools::check_win_devel()`

## R CMD check results

There were no ERRORs or WARNINGs.
There was one NOTE:

> The Description field should not start with the package name,
  'This package' or similar.

False positive: the adjective 'Rogue' happens to match the name of the package.
