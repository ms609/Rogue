#' @param fullSeq Logical specifying whether to list all taxa (`TRUE`), or
#' only those that improve information content when all are dropped (`FALSE`).
#' @inheritParams TipInstability
#' @describeIn RogueTaxa Shortcut to 'fast' heuristic, with option to return
#' evaluation of all taxa using `fullSeq = TRUE`.
#' @examples
#'
#' QuickRogue(trees, fullSeq = TRUE)
#' @importFrom ape consensus
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom TreeDist ConsensusInfo
#' @importFrom fastmatch %fin%
#' @importFrom TreeTools NTip SplitFrequency
#' @export
QuickRogue <- function (trees,
                        info = 'phylogenetic',
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
                        IC = ConsensusInfo(c(trees), info = info),
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
  nTree <- length(trees)

  tr <- trees
  neverDrop <- .NeverDrop(neverDrop, trees[[1]]$tip.label)
  nKeep <- length(neverDrop)
  candidates <- character(nTip - 2L - nKeep)
  score <- double(nTip - 2)
  score[1] <- ConsensusInfo(trees, info = info, check.tips = FALSE)
  nDrops <- nTip - 3L - nKeep
  cli_progress_bar("Dropping leaves", total = nDrops * (nDrops + 1L) / 2 )
  for (i in 1 + seq_len(nDrops)) {
    cli_progress_update(nDrops - (i - 1),
                        status = paste0("Leaf ", i - 1, "/", nDrops))
    tipScores <- TipInstability(tr, log = log, average = average,
                                deviation = deviation,
                                checkTips = FALSE)
    tipScores[tr[[1]]$tip.label %fin% neverDrop] <- -Inf
    candidate <- which.max(tipScores)
    if (length(candidate)) {
      candidates[i] <- names(candidate)
    }
    tr <- lapply(tr, DropTip, candidate, preorder = FALSE)
    score[i] <- ConsensusInfo(tr, info = info, check.tips = FALSE)
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
             taxNum = c(NA_character_, fmatch(dropped, trees[[1]]$tip.label)),
             taxon = c(NA_character_, dropped),
             rawImprovement = c(NA_real_, score[-1] - score[-length(score)]),
             IC = score,
             stringsAsFactors = FALSE)
}

#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' cli_alert_success
#' @importFrom fastmatch fmatch
#' @importFrom TreeDist ConsensusInfo
#' @importFrom TreeTools DropTip SplitFrequency Preorder RenumberTips
#' @importFrom utils combn
Roguehalla <- function (trees, dropsetSize = 1, info = 'phylogenetic',
                        neverDrop) {
  if (!inherits(trees, 'multiPhylo')) {
    if (inherits(trees, 'phylo')) {
      return(data.frame(num = 0,
                        taxNum = NA_character_,
                        taxon = NA_character_,
                        rawImprovement = NA_real_,
                        IC = ConsensusInfo(c(trees), info = info),
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
  nTree <- length(trees)
  majority <- 0.5 + sqrt(.Machine$double.eps)

  startTip <- NTip(trees[[1]])
  neverDrop <- .NeverDrop(neverDrop, trees[[1]]$tip.label)
  best <- ConsensusInfo(trees, info = info, check.tips = FALSE)

  .Drop <- function (n) {
    cli_progress_bar(paste0("Dropset size ", n))
    keepN <- fmatch(neverDrop, trees[[1]]$tip.label)
    nTip <- NTip(trees[[1]])
    nKept <- nTip - length(keepN)
    drops <- matrix(setdiff(seq_len(nTip), keepN)[combn(nKept, n)], nrow = n)

    cli_progress_update(set = 0, total = ncol(drops))
    candidates <- apply(drops, 2, function (drop) {
      cli_progress_update(1, .envir = parent.frame(2), status = paste0(
        "Drop ", startTip - NTip(trees[[1]]), " leaves = ",
        signif(best), " bits."))
      dropForest <- lapply(trees, DropTip, drop, preorder = FALSE)
      ConsensusInfo(dropForest, info = info, check.tips = FALSE)
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
        trees <- lapply(trees, DropTip, thisDrop, preorder = FALSE)
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
