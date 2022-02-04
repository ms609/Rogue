#' @param fullSeq Logical specifying whether to list all taxa (`TRUE`), or
#' only those that improve information content when all are dropped (`FALSE`).
#' @param p Proportion of trees that must contain a split before it is included
#' in the consensus under consideration.  0.5, the default, corresponds to a
#' majority rule tree; 1.0 will maximize the information content of the
#' strict consensus.
#' @inheritParams TipInstability
#' @describeIn RogueTaxa Shortcut to 'fast' heuristic, with option to return
#' evaluation of all taxa using `fullSeq = TRUE`.
#' @examples
#'
#' QuickRogue(trees, fullSeq = TRUE)
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom TreeDist ConsensusInfo SplitwiseInfo ClusteringInfo
#' @importFrom fastmatch %fin%
#' @importFrom TreeTools NTip SplitFrequency PectinateTree DropTipPhylo
#' @export
QuickRogue <- function(trees,
                        info = 'phylogenetic',
                        p = 0.5,
                        log = TRUE, average = 'median', deviation = 'mad',
                        neverDrop, fullSeq = FALSE) {
  if (!is.na(pmatch(tolower(info), 'spic'))) {
    info <- 'phylogenetic'
  } else if (!is.na(pmatch(tolower(info), 'scic'))) {
    info <- 'clustering'
  }
  if (!inherits(trees, 'multiPhylo')) {
    if (inherits(trees, 'phylo')) {
      return(data.frame(num = 0,
                        taxNum = NA_character_,
                        taxon = NA_character_,
                        rawImprovement = NA_real_,
                        IC = ConsensusInfo(c(trees), info = info, p = p),
                        stringsAsFactors = FALSE))
    }
    if (!is.list(trees)) {
      stop("`trees` must be a list of `phylo` objects")
    }
    trees <- structure(trees, class = 'multiPhylo')
  }
  trees <- lapply(trees, RenumberTips, trees[[1]])
  trees <- lapply(trees, Preorder)
  nTip <- NTip(trees[[1]])

  tr <- trees
  neverDrop <- .NeverDrop(neverDrop, trees[[1]]$tip.label)
  nKeep <- length(neverDrop)
  candidates <- character(nTip - 2L - nKeep)
  score <- double(nTip - 2)
  score[1] <- ConsensusInfo(trees, info = info, p = p, check.tips = FALSE)
  nDrops <- nTip - 3L - nKeep
  TotalInfo <- switch(pmatch(info, c('phylogenetic', 'clustering')),
                      SplitwiseInfo, ClusteringInfo)

  cli_progress_bar("Drop leaf", total = nDrops * (nDrops + 1L) / 2,
                   .auto_close = FALSE)
  for (i in 1L + seq_len(nDrops)) {
    bestPossibleNext <- TotalInfo(PectinateTree(nTip - i))
    bestYet <- max(score, na.rm = TRUE)
    if (bestPossibleNext < bestYet) {
      # message("Broken out: can't attain ", signif(score[i]), " bits with ",
      #         nTip - i, " leaves.")
      break
    }
    bitStat <- if (bestPossibleNext < 1000) {
      paste0(round(bestYet), '/', round(bestPossibleNext), ' bits')
    } else if (bestPossibleNext < 1000000) {
      paste0(round(bestYet / 1000), '/', round(bestPossibleNext / 1000), ' kb')
    } else {
      paste0(round(bestYet / 1e6), '/', round(bestPossibleNext / 1e6), ' Mb')
    }
    cli_progress_update(nDrops - (i - 1),
                        status = paste0("Leaf ", i - 1, '; ', bitStat))
    tipScores <- TipInstability(tr, log = log, average = average,
                                deviation = deviation,
                                checkTips = FALSE)
    tipScores[tr[[1]]$tip.label %fin% neverDrop] <- -Inf
    candidate <- which.max(tipScores)
    if (length(candidate)) {
      candidates[i] <- names(candidate)
    }
    tr <- lapply(tr, DropTipPhylo, candidate, preorder = FALSE, check = FALSE)
    score[i] <- ConsensusInfo(tr, info = info, p = p, check.tips = FALSE)
  }
  cli_progress_done()

  bestPos <- which.max(score)
  bestScore <- score[bestPos]
  pointer <- bestPos - 1L
  needsRecalc <- logical(length(candidates))

  cli_progress_bar("Restore leaf", total = bestPos - 2L)
  while (pointer > 1L) {
    tryScore <- ConsensusInfo(
      lapply(trees, DropTipPhylo, candidates[seq_len(bestPos)[-c(1, pointer)]],
             preorder = FALSE, check = FALSE),
      info = info, p = p, check.tips = FALSE)
    if (tryScore > bestScore) {
      candidates[1:bestPos] <- candidates[c((1:bestPos)[-pointer], pointer)]
      bestScore <- tryScore
      bestPos <- bestPos - 1L
      score[bestPos] <- tryScore
      needsRecalc[bestPos - seq_len(bestPos - pointer)] <- TRUE
    }
    cli_progress_update(1)
    pointer <- pointer - 1L
  }
  for (i in which(needsRecalc)) {
    score[i] <- ConsensusInfo(lapply(trees, DropTipPhylo,
                                     preorder = FALSE, check = FALSE,
                                     candidates[seq_len(i)[-1]]),
                              info = info, p = p, check.tips = FALSE)
  }
  cli_progress_done()

  dropped <- if (fullSeq) {
    union(candidates, trees[[1]]$tip.label)[-1]
  } else {
    candidates[seq_len(which.max(score))[-1]]
  }
  score <- score[seq_len(length(dropped) + 1L)]
  score[is.na(score)] <- 0

  # Return:
  data.frame(num = seq_along(score) - 1L,
             taxNum = c(NA_character_, match(dropped, trees[[1]]$tip.label)),
             taxon = c(NA_character_, dropped),
             rawImprovement = c(NA_real_, score[-1] - score[-length(score)]),
             IC = score,
             stringsAsFactors = FALSE)
}

#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' cli_alert_success
#' @importFrom fastmatch fmatch
#' @importFrom TreeDist ConsensusInfo
#' @importFrom TreeTools DropTipPhylo SplitFrequency Preorder RenumberTips
#' @importFrom utils combn
Roguehalla <- function(trees, dropsetSize = 1, info = 'phylogenetic',
                        p = 0.5, neverDrop) {
  if (!inherits(trees, 'multiPhylo')) {
    if (inherits(trees, 'phylo')) {
      return(data.frame(num = 0,
                        taxNum = NA_character_,
                        taxon = NA_character_,
                        rawImprovement = NA_real_,
                        IC = ConsensusInfo(c(trees), info = info, p = p),
                        stringsAsFactors = FALSE))
    }
    if (!is.list(trees)) {
      stop("`trees` must be a list of `phylo` objects")
    }
    trees <- structure(trees, class = 'multiPhylo')
  }
  trees <- lapply(trees, RenumberTips, trees[[1]])
  trees <- lapply(trees, Preorder)
  startTrees <- trees
  labels <- startTrees[[1]]$tip.label

  startTip <- NTip(trees[[1]])
  neverDrop <- .NeverDrop(neverDrop, trees[[1]]$tip.label)
  best <- ConsensusInfo(trees, info = info, p = p, check.tips = FALSE)

  .Drop <- function(n) {
    cli_progress_bar(paste0("Dropset size ", n))
    keepN <- fmatch(neverDrop, trees[[1]]$tip.label)
    nTip <- NTip(trees[[1]])
    nKept <- nTip - length(keepN)
    drops <- matrix(setdiff(seq_len(nTip), keepN)[combn(nKept, n)], nrow = n)

    cli_progress_update(set = 0, total = ncol(drops))
    candidates <- apply(drops, 2, function(drop) {
      cli_progress_update(1, .envir = parent.frame(2), status = paste0(
        "Drop ", startTip - NTip(trees[[1]]), " leaves = ",
        signif(best), " bits."))
      dropForest <- lapply(trees, DropTipPhylo, drop,
                           check = FALSE, preorder = FALSE)
      ConsensusInfo(dropForest, info = info, p = p, check.tips = FALSE)
    })
    cli_progress_done()
    if (max(candidates) > best) {
      list(info = max(candidates),
           drop = trees[[1]]$tip.label[drops[, which.max(candidates)]])
    } else {
      NULL
    }
  }

  dropSeq <- integer(0)
  taxSeq <- character(0)
  dropInf <- best
  repeat {
    improved <- FALSE
    for (i in seq_len(dropsetSize)) {
      dropped <- .Drop(i)
      if (!is.null(dropped)) {
        improved <- TRUE
        best <- dropped$info
        thisDrop <- dropped[['drop']]
        dropSeq <- c(dropSeq, paste0(thisDrop, collapse = ','))
        taxSeq <- c(taxSeq, paste0(fmatch(thisDrop, labels), collapse = ','))
        dropInf <- c(dropInf, best)
        trees <- lapply(trees, DropTipPhylo, thisDrop, preorder = FALSE,
                        check = FALSE)
        break
      }
    }
    if (!improved) {
      cli_alert_success(paste0("{Sys.time()}: Dropped ",
                               startTip - NTip(trees[[1]]),
                               " leaves, rendering {signif(best)} bits."))
      break
    }
  }

  # Return:
  data.frame(num = c(0, seq_along(dropSeq)),
             taxNum = c(NA_character_, taxSeq),
             taxon = c(NA_character_, dropSeq),
             rawImprovement = c(NA, dropInf[-1] - dropInf[-length(dropInf)]),
             IC = dropInf,
             stringsAsFactors = FALSE)
}
