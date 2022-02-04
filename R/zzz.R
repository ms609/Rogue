.onUnload <- function(libpath) {
  library.dynam.unload("Rogue", libpath)
}

# Additional tests:
#
# spell_check()
#
# run_examples()
#
# devtools::check_win_devel(quiet = TRUE); rhub::check_for_cran()
# Check valgrind results on Github Actions
# revdepcheck::revdep_check()
#
# codemetar::write_codemeta()
