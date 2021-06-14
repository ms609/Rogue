#' Call RogueNaRok
#' @useDynLib RogueTaxa, .registration = TRUE
#' @template MRS
#' @export
RogueNaRok <- function (treeFile = "",
                        bootTrees = "", 
                        computeSupport = FALSE,
                        run_id = "",
                        maxDropsetSize = 0,
                        excludeFile = "",
                        workdir = "",
                        labelPenalty = "",
                        mreOptimization = FALSE,
                        threshold = 50
                        ) {
  .Call("RogueNaRok",
        R_bootTrees = as.character(bootTrees),
        R_computeSupport = as.logical(computeSupport),
        R_run_id = as.character(run_id),
        R_treeFile = as.character(treeFile),
        R_maxDropsetSize = as.double(maxDropsetSize),
        R_excludeFile = as.character(excludeFile),
        R_workdir = as.character(workdir),
        R_labelPenalty = as.double(labelPenalty),
        R_mreOptimization = as.logical(mreOptimization),
        R_threshold = as.double(threshold))
}