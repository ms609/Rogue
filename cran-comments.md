## Test environments

* Local PC:
  - Windows 10, R 4.3.1

* [GitHub Actions](https://github.com/ms609/Rogue/actions)
  - Ubuntu 20.04
    - R 4.1
    - R release (tests, examples & vignettes run with valgrind)
    - R devel
  - Mac OS X 12.6.8, R release
  - Microsoft Windows Server 2022 10.0.20348, R release
  
* R-hub, with `rhub::check_for_cran()` and `devtools::check_win_devel()`

## R CMD check results

There were no ERRORs or WARNINGs.
There was one NOTE:

> The Description field should not start with the package name,
  'This package' or similar.

False positive: the word 'Rogue' is used here as an adjective, but happens to
match the name of the package.
