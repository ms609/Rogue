#' Drop rogue taxa to generate a more informative consensus
#'
#' `RogueTaxa()` finds wildcard leaves whose removal increases the resolution
#' or branch support values of a consensus tree, using the relative
#' bipartition, shared phylogenetic, or mutual clustering concepts of
#' information.
#'
#' "Rogue" or (loosely) "wildcard" taxa \insertCite{Nixon1992}{Rogue} are
#' leaves whose position in a tree is poorly constrained, typically because
#' much of the phylogenetic data associated with the taxon is either missing or
#' in conflict with other data \insertCite{Kearney2002}{Rogue}.
#'
#' These functions use heuristic methods to identify rogue taxa whose removal
#' improves the information content of a consensus tree, by the definitions
#' of information discussed below.
#'
#' @section Information criteria:
#' The splitwise phylogenetic information content measure produces the best
#' results \insertCite{SmithCons}{Rogue}.
#' It uses the splitwise information content as a shortcut, which involves
#' double counting of some information (which may or may not be desirable).
#' The same holds for the mutual clustering information measure; this measure
#' is less obviously suited to the detection of rogues.
#'
#' The "relative bipartition information criterion" (\acronym{RBIC}) is
#' the sum of all support values divided by the maximum possible support in a
#' fully bifurcating tree with the initial set of taxa.
#' The relative bipartition information content approach employs the
#' 'RogueNaRok' implementation \insertCite{Aberer2013}{Rogue}, which can handle
#' large trees relatively quickly.
#' The \acronym{RBIC} is is not strictly a measure of information and can
#' produce undesirable results \insertCite{Wilkinson2017}{Rogue}.
#'
#' `C_RogueNaRok()` directly interfaces the 'RogueNaRok' C implementation,
#' with no input checking; be aware that invalid input will cause undefined
#' behaviour and is likely to crash R.
#'
#' @param trees List of trees to analyse.
#' @param info Concept of information to employ; see details.
#' @param neverDrop Tip labels that should not be dropped from the consensus.
#' @param verbose Logical specifying whether to display output from RogueNaRok.
#' If `FALSE`, output will be included as an attribute of the return value.
#' @param return If `taxa`, returns the leaves identified as rogues; if `tree`,
#' return a consensus tree omitting rogue taxa.
#'
#' @return `RogueTaxa()` returns a `data.frame`. Each row after the first,
#' which describes the starting tree, describes a dropset operation.
#' Columns describe:
#' - `num`: Sequential index of the drop operation
#' - `taxNum`: Numeric identifier of the dropped leaves
#' - `taxon`: Text identifier of dropped leaves
#' - `rawImprovement`: Improvement in score obtained by this operation
#' - `IC`: Information content of tree after dropping all leaves so far,
#' by the measure indicated by `info`.
#'
#' @examples
#' library("TreeTools", warn.conflicts = FALSE)
#'
#' trees <- list(read.tree(text = ("(a, (b, (c, (d, (e, (X1, X2))))));")),
#'               read.tree(text = ("((a, (X1, X2)), (b, (c, (d, e))));")))
#' RogueTaxa(trees, dropsetSize = 2)
#'
#' trees <- list(
#'      read.tree(text = '((a, y), (b, (c, (z, ((d, e), (f, (g, x)))))));'),
#'      read.tree(text = '(a, (b, (c, (z, (((d, y), e), (f, (g, x)))))));'),
#'      read.tree(text = '(a, (b, ((c, z), ((d, (e, y)), ((f, x), g)))));'),
#'      read.tree(text = '(a, (b, ((c, z), ((d, (e, x)), (f, (g, y))))));'),
#'      read.tree(text = '(a, ((b, x), ((c, z), ((d, e), (f, (g, y))))));')
#'      )
#' cons <- consensus(trees, p = 0.5)
#' plot(cons)
#' LabelSplits(cons, SplitFrequency(cons, trees) / length(trees))
#' reduced <- RogueTaxa(trees, info = 'phylogenetic', ret = 'tree')
#' plot(reduced)
#' LabelSplits(reduced, SplitFrequency(reduced, trees) / length(trees))
#' @author [Martin R. Smith](https://smithlabdurham.github.io/)
#' (<martin.smith@durham.ac.uk>), linking to
#' [RogueNaRok](https://github.com/aberer/RogueNaRok/)
#' C library by Andre Aberer (<andre.aberer at googlemail.com>)
#'
#' @importFrom ape write.tree
#' @importFrom TreeTools Consensus ConsensusWithout RenumberTips
#' @importFrom TreeDist SplitwiseInfo ClusteringInfo ConsensusInfo
#' @importFrom utils capture.output read.table
#' @references \insertAllCited{}
#' @export
RogueTaxa <- function (trees,
                       info = c('spic', 'scic', 'fspic', 'fscic', 'rbic'),
                       return = c('taxa', 'tree'),
                       bestTree = NULL,
                       computeSupport = TRUE,
                       dropsetSize = 1,
                       neverDrop = character(0),
                       labelPenalty = 0,
                       mreOptimization = FALSE,
                       threshold = 50,
                       verbose = FALSE) {
  p <- threshold / 100

  # Check format of `trees`
  if (!inherits(trees, 'multiPhylo')) {
    if (inherits(trees, 'phylo')) {
      return(data.frame(num = 0,
                        taxNum = NA_character_,
                        taxon = NA_character_,
                        rawImprovement = NA_real_,
                        IC = ConsensusInfo(c(trees), info = info[1], p = p),
                        stringsAsFactors = FALSE))
    }
    trees <- structure(trees, class = 'multiPhylo')
  }

  labels <- attr(trees, 'tip.label')
  if (is.null(labels)) {
    labels <- lapply(trees, `[[`, 'tip.label')
  }
  if (is.list(labels)) {
    if (length(unique(vapply(labels, length, 1))) > 1) {
      stop("All trees must bear the same number of labels.");
    }
    leaves <- labels[[1]]
    lapply(labels[-1], function (these) {
      if (length(setdiff(leaves, these))) {
        stop("All trees must bear the same labels.\n  Found tree lacking ",
             paste0(setdiff(leaves, these), collapse = ', '), '.')
      }
    })
  } else {
    leaves <- labels
  }
  if (any(duplicated(leaves))) {
    stop("All leaves must bear unique labels.")
  }
  trees <- lapply(trees, RenumberTips, trees[[1]])
  trees <- structure(lapply(trees, Preorder), class = 'multiPhylo')

  # Select and apply method
  info <- tolower(info[1])
  if (!is.na(pmatch(info, 'phylogenetic'))) {
    info <- 2L
  } else if (!is.na(pmatch(info, 'clustering'))) {
    info <- 3L
  } else {
    info <- pmatch(info, c('rbic', 'spic', 'scic', 'fspic', 'fscic'))
    if (is.na(info)) {
      stop("`info` must be 'rbic', 'spic', 'fspic', 'scic' or 'fscic'")
    }
  }

  result <- switch(info,
                   .RogueNaRok(trees, bestTree = bestTree,
                               computeSupport = computeSupport,
                               dropsetSize = dropsetSize, neverDrop = neverDrop,
                               labelPenalty = labelPenalty,
                               mreOptimization = mreOptimization,
                               threshold = threshold, verbose = verbose),
                   Roguehalla(trees, dropsetSize = dropsetSize, info = 'phylo',
                              p = p, neverDrop = neverDrop),
                   Roguehalla(trees, dropsetSize = dropsetSize, info = 'clust',
                              p = p, neverDrop = neverDrop),
                   QuickRogue(trees, info = 'phylo', p = p,
                              neverDrop = neverDrop),
                   QuickRogue(trees, info = 'clust', p = p,
                              neverDrop = neverDrop)
  )

  # Format return value
  returnTree <- 1L == pmatch(tolower(return),
                             c('tree', 'taxa', 'tips', 'leaves'))[1]
  if (is.na(returnTree)) {
    stop('`return` must be "tree" or "tips".')
  }

  # Return:
  if (returnTree) {
    drops <- unlist(strsplit(result[-1, 'taxon'], ','))
    if (is.null(drops)) {
      Consensus(trees, p = p)
    } else {
      ConsensusWithout(trees, drops, p = p)
    }
  } else {
    result
  }
}
