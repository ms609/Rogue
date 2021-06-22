.onUnload <- function (libpath) {
  library.dynam.unload("Rogue", libpath)
}

# Additional tests:
#
# spell_check()
# pkgdown::build_reference_index()
#
# run_examples()
# build_vignettes()
#
# devtools::check_win_devel(quiet = TRUE); rhub::check_for_cran()
# Check valgrind results on Github Actions
# revdepcheck::revdep_check()
#
# codemetar::write_codemeta()
#
# tools::resaveRdaFiles('R', compress='auto') - is default bzip2 the optimal?
# tools::checkRdaFiles('R') - set optimal compression in `data-raw`
