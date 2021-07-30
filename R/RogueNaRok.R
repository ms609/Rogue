.RogueNaRok <- function (trees,
                         bestTree = NULL,
                         computeSupport = TRUE,
                         dropsetSize = 1,
                         neverDrop = character(0),
                         labelPenalty = 0,
                         mreOptimization = FALSE,
                         threshold = 50,
                         verbose = FALSE) {
  wd <- tempdir()

  bootTrees <- tempfile(tmpdir = wd)
  write.tree(trees, file = bootTrees)
  on.exit(unlink(bootTrees))


  if (inherits(bestTree, 'phylo')) {
    treeFile <- tempfile(tmpdir = wd)
    write.tree(bestTree, treeFile)
    on.exit(unlink(treeFile))

  } else {
    treeFile <- ""
  }
  if (length(neverDrop)) {
    excludeFile <- tempfile(tmpdir = wd)
    write(neverDrop, excludeFile)
    on.exit(unlink(excludeFile))
  } else {
    excludeFile <- ""
  }
  RunRogueNaRok <- function ()
    C_RogueNaRok(bootTrees = bootTrees,
                 runId = "tmp",
                 treeFile = treeFile,
                 computeSupport = computeSupport,
                 dropsetSize = dropsetSize,
                 excludeFile = excludeFile,
                 workDir = wd,
                 labelPenalty = labelPenalty,
                 mreOptimization = mreOptimization,
                 threshold = threshold)
  if (verbose) {
    RunRogueNaRok()                                                             # nocov
  } else {
    output <- capture.output(RunRogueNaRok())
  }


  rogueFile <- paste0(wd, '/RogueNaRok_droppedRogues.tmp')
  if (!file.exists(rogueFile)) {
    stop("RogueNaRok did not produce output at ", rogueFile)
  }
  droppedRogues <- read.table(rogueFile, header = TRUE,
                              colClasses = c('integer', 'character',
                                             'character', 'numeric', 'numeric'))

  unlink(rogueFile)
  unlink(paste0(wd, '/RogueNaRok_info.tmp'))
  if (verbose) {
    droppedRogues                                                               # nocov
  } else {
    structure(droppedRogues, log = output)
  }
}

#' @param bootTrees A collection of bootstrap trees.
#' @param runId An identifier for this run, appended to output files.
#' @param treeFile,bestTree If a single best-known tree (such as an ML or MP tree)
#' is provided, RogueNaRok optimizes the bootstrap support in this
#' best-known tree (still drawn from the bootstrap trees).
#' The `threshold` parameter is ignored.
#' @param excludeFile Taxa in this file (one taxon per line) will not be
#' considered for pruning.
#' @param threshold,mreOptimization A threshold or mode for the consensus tree
#' that is optimized. Specify a value between 50 (majority rule consensus) and
#' 100 (strict consensus), or set `mreOptimization = TRUE`
#' for the extended majority rule consensus.
#' Note that rogue taxa identified with respect to different thresholds can
#' vary substantially. DEFAULT: MR consensus.
#' @param computeSupport Logical: Instead of trying to maximize the support
#' in the consensus tree, the RogueNaRok will try to maximize the number of
#' bipartitions in the final tree by pruning taxa.
#' @param labelPenalty A weight factor to penalize for dropset size.
#' `labelPenalty = 1` is Pattengale's criterion.
#' The higher the value, the more conservative the algorithm is in pruning taxa.
#' DEFAULT: 0.0 (=RBIC)
#' `info = 'rbic'` only; other measures implicitly penalize dropset size.
#' @param dropsetSize Maximum size of dropset per iteration.
#' If `dropsetSize == n`, then RogueNaRok will test in each iteration which
#' tuple of n taxa increases optimality criterion the most and prunes
#' taxa accordingly.
#' This improves the result, but run times will increase at least linearly.
#' @param workDir A working directory where output files are created.
#' @return `C_RogueNaRok()` returns `0` if successful; `-1` on error.
#' @useDynLib Rogue, .registration = TRUE
#' @rdname RogueTaxa
#' @examples
#'
#' bootTrees <- system.file('example/150.bs', package = 'Rogue')
#' tmpDir <- tempdir()
#' C_RogueNaRok(bootTrees, workDir = tmpDir)
#'
#' # Results have been written to our temporary directory
#' oldwd <- setwd(tmpDir)
#' read.table('RogueNaRok_droppedRogues.tmp', header = TRUE)
#'
#' # Delete temporary files
#' file.remove('RogueNaRok_droppedRogues.tmp')
#' file.remove('RogueNaRok_info.tmp')
#'
#' setwd(oldwd)
#' @export
C_RogueNaRok <- function (bootTrees = "",
                          runId = "tmp",
                          treeFile = "",
                          computeSupport = TRUE,
                          dropsetSize = 1,
                          excludeFile = "",
                          workDir = "",
                          labelPenalty = 0,
                          mreOptimization = FALSE,
                          threshold = 50
                          ) {
  .Call("RogueNaRok",
        R_bootTrees = as.character(bootTrees),
        R_run_id = as.character(runId),
        R_treeFile = as.character(treeFile),
        R_computeSupport = as.logical(computeSupport),
        R_maxDropsetSize = as.integer(dropsetSize),
        R_excludeFile = as.character(excludeFile),
        R_workdir = as.character(workDir),
        R_labelPenalty = as.double(labelPenalty),
        R_mreOptimization = as.logical(mreOptimization),
        R_threshold = as.double(threshold))
}
