## Test environments

* Local PC:
  - Windows 10, R 4.1.0-patched

* [GitHub Actions](https://github.com/ms609/Roguer/actions)
  - Ubuntu 20.04
    - R 3.6.3
    - R release (tests, examples & vignettes run with valgrind)
    - R devel
  - Mac OS X 10.15.7, R release
  - Microsoft Windows Server 2019 10.0.17763, R release
  
* R-hub, with `rhub::check_for_cran()` and `devtools::check_win_devel()`

## R CMD check results

There were no ERRORs or WARNINGs.
There was one NOTE:

> Possibly mis-spelled words in DESCRIPTION:
>   RogueNaRok (14:27)

These spellings have been verified.
