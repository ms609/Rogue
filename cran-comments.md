## Test environments

* Local PC:
  - Windows 10, R 4.1.1

* [GitHub Actions](https://github.com/ms609/Rogue/actions)
  - Ubuntu 20.04
    - R 3.6.3
    - R release (tests, examples & vignettes run with valgrind)
    - R devel
  - Mac OS X 10.15.7, R release
  - Microsoft Windows Server 2019 10.0.17763, R release
  
* R-hub, with `rhub::check_for_cran()` and `devtools::check_win_devel()`

## R CMD check results

There were no ERRORs or WARNINGs.
There were two NOTEs:

> Possibly mis-spelled words in DESCRIPTION:
>   RogueNaRok (14:27)

These spellings are correct.

> The Description field should not start with the package name,
  'This package' or similar.

False positive: the word 'Rogue' is used here as an adjective, but happens to
match the name of the package.
