
#' @export
RogueNaRok <- function (trees, bestTree = NULL, 
                        computeSupport = TRUE,
                        maxDropsetSize = 1,
                        neverDrop = character(0),
                        labelPenalty = 0,
                        mreOptimization = FALSE,
                        threshold = 50) {
  wd <- tempdir()
  bootTrees <- tmpfile(write.tree(trees), tmpdir = wd)
  treeFile <- if (inherits(bestTree, 'phylo')) {
    tmpfile(write.tree(bestTree), tmpdir = wd)
  } else {
    ""
  }
  excludeFile <- if (length(neverDrop)) {
    tmpfile(neverDrop, tmpdir = wd)
  } else {
    ""
  }
  C_RogueNaRok(bootTrees = bootTrees, 
               run_id = "rnr.tmp",
               treeFile = treeFile,
               computeSupport = computeSupport,
               maxDropsetSize = maxDropsetSize,
               excludeFile = excludeFile,
               workdir = wd,
               labelPenalty = labelPenalty,
               mreOptimization = mreOptimization,
               threshold = threshold)
  rogueFile <- paste0(wd, '\\RougeNaRokR_droppedRogues.rnr.tmp')
  droppedRogues <- read.table(rogueFile, header = TRUE)
  file.remove(rogueFile)
  file.remove(paste0(wd, '\\RougeNaRokR_info.rnr.tmp'))
  droppedRogues
}

#' Call RogueNaRok
#'
#' Implements the RogueNaRok algorithm for rogue taxon identification.
#'  
#' @param bootTrees A collection of bootstrap trees.
#' @param runId An identifier for this run.
#' @param treeFile If a single best-known tree (such as an ML or MP tree) 
#' is provided, RogueNaRok optimizes the bootstrap support in this
#' best-known tree (still drawn from the bootstrap trees).
#' The `threshold` parameter is ignored.
#' @param excludeFile Taxa in this file (one taxon per line) will not be
#' considered for pruning.
#' @param threshold A threshold or mode for the consensus tree that is
#' optimized. Specify a value between 50 (majority rule consensus) and
#' 100 (strict consensus) or MR (for the extended majority rule
#' consensus). Note that rogue taxa identified with respect to different
#' thresholds can vary substantially. DEFAULT: MR consensus
#' @param computeSupport Logical: Instead of trying to maximize the support
#' in the consensus tree, the RogueNaRok will try to maximize the number of 
#' bipartitions in the final tree by pruning taxa.
#' @param labelPenalty A weight factor to penalize for dropset size.
#' `labelPenalty = 1`  is Pattengale's criterion. 
#' The higher the value, the more conservative the algorithm is in pruning taxa.
#' DEFAULT: 0.0 (=RBIC)
#' @param dropsetSize Maximum size of dropset per iteration. 
#' If `dropsetSize == n`, then RogueNaRok will test in each iteration which
#' tuple of n taxa increases optimality criterion the most and prunes
#' taxa accordingly. 
#' This improves the result, but runtimes will increase at least linearly.
#' @param workDir A working directory where output files are created.
#' @return `C_RogueNaRok()` returns `0` if successful; `-1` on error.
#' @useDynLib RogueTaxa, .registration = TRUE
#' @template MRS
#' @rdname RogueNaRok
#' @export
C_RogueNaRok <- function (bootTrees = "", 
                          run_id = "tmp",
                          treeFile = "",
                          computeSupport = TRUE,
                          maxDropsetSize = 1,
                          excludeFile = "",
                          workdir = "",
                          labelPenalty = 0,
                          mreOptimization = FALSE,
                          threshold = 50
                          ) {
  .Call("RogueNaRok",
        R_bootTrees = as.character(bootTrees),
        R_run_id = as.character(run_id),
        R_treeFile = as.character(treeFile),
        R_computeSupport = as.logical(computeSupport),
        R_maxDropsetSize = as.double(maxDropsetSize),
        R_excludeFile = as.character(excludeFile),
        R_workdir = as.character(workdir),
        R_labelPenalty = as.double(labelPenalty),
        R_mreOptimization = as.logical(mreOptimization),
        R_threshold = as.double(threshold))
}
