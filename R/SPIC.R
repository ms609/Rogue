#' @importFrom stats cophenetic
Cophenetic <- function (x) {
        if (is.null(x$edge.length)) {
                x$edge.length <- rep_len(1, dim(x$edge)[1])
        }
        ret <- cophenetic(x)
        ret[ret < sqrt(.Machine$double.eps)] <- 0
        ret
}

# @examples
# library("TreeTools", quietly = TRUE)
# trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')
# plot(consensus(trees))
# instab <- TipInstability(trees)
# plot(ConsensusWithout(trees, names(instab[instab > 0.2])))
# @template MRS
#' @importFrom stats cmdscale mad
TipInstability <- function (trees) {
  dists <- .TipDistances(trees)

  means <- rowMeans(dists, dims = 2)
  devs <- apply(dists, 1:2, function(x) mad(x))
  diag(devs) <- NA
  relDevs <- devs / mean(means[lower.tri(means)])
  rowMeans(relDevs, na.rm = TRUE)
}

.TipDistances <- function (trees) {
  nTip <- NTip(trees)
  if (length(unique(nTip)) > 1) {
    stop("Trees must have same number of leaves")
  }
  nTip <- nTip[1]
  trees[-1] <- lapply(trees[-1], RenumberTips, trees[[1]])
  dists <- vapply(trees, Cophenetic, matrix(0, nTip, nTip))
}

#' @importFrom grDevices hcl
#' @importFrom stats setNames
.TipCols <- function (trees, luminence = 50) {
  dists <- .TipDistances(trees)

  means <- rowMeans(dists, dims = 2)
  devs <- apply(dists, 1:2, function(x) mad(x))
  diag(devs) <- NA
  relDevs <- devs / mean(means[lower.tri(means)])

  pc <- cmdscale(means, k = 1)
  pc <- pc - min(pc)
  pc <- pc * 340 / max(pc)

  setNames(hcl(h =  pc, c = 100 * (1 - rowMeans(relDevs, na.rm = TRUE)),
               l = luminence),
           TipLabels(trees[[1]]))

}

#' Detect rogue taxa using phylogenetic information distance
#'
#' Calculate the volatility of each tip: namely, the impact on the mean
#' phylogenetic information distance between trees when that tip is removed.
#'
#' @inheritParams RogueTaxa
#' @return `TipVolatility()` returns a named vector listing the volatility
#' index calculated for each leaf.
#' @references
#' \insertRef{Aberer2013}{Rogue}
#' \insertRef{Wilkinson2017}{Rogue}
#' @examples
#' library("TreeTools", quietly = TRUE)
#' trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')
#' trees[] <- lapply(trees, AddTip, 'Rogue', 'Rogue2')
#'
#' TipInformation(trees)
#' sb <- TipInformation(trees)
#' col <- hcl.colors(ceiling(max(sb) *1.13), 'inferno')[ceiling(sb)]
#' plot(consensus(trees), tip.color = col)
#' plot(ConsensusWithout(trees, names(sb[sb == max(sb)])))
#' @importFrom TreeDist PhylogeneticInfoDistance
#' @importFrom TreeTools CladisticInfo
#' @export
TipVolatility <- function (trees) {
  tips <- trees[[1]]$tip.label
  startInfo <- mean(CladisticInfo(trees))
  info <- vapply(tips, function (drop) {
          tr <- lapply(trees, DropTip, drop)
          c(meanInfo = mean(CladisticInfo(tr)),
            meanDist = mean(PhylogeneticInfoDistance(tr, normalize = TRUE)))
  }, double(2))
  mean(PhylogeneticInfoDistance(trees, normalize = TRUE)) - info['meanDist', ]
}

#' @param fullSeq Logical specifying whether to list: `TRUE`, all taxa, or
#' `FALSE`, only those that improve information content when all are dropped.
#' @describeIn RogueTaxa Shortcut to 'fast' heuristic, with option to return
#' evaluation of all taxa using `fullSeq = TRUE`.
#' @examples
#'
#' QuickRogue(trees, fullSeq = TRUE)
#' @importFrom ape consensus
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom TreeDist ConsensusInfo
#' @importFrom TreeTools NTip SplitFrequency
#' @export
QuickRogue <- function (trees, info = 'clustering', fullSeq = FALSE) {
  if (!is.na(pmatch(tolower(info), 'spic'))) {
    info <- 'phylogenetic'
  } else if (!is.na(pmatch(tolower(info), 'scic'))) {
    info <- 'clustering'
  }
  if (!inherits(trees, 'multiPhylo')) {
    if (inherits(trees, 'phylo')) {
      return (trees)
    }
    if (!is.list(trees)) {
      stop("`trees` must be a list of `phylo` objects")
    }
    trees <- structure(trees, class = 'multiPhylo')
  }
  lastMessage <- Sys.time()
  trees <- lapply(trees, RenumberTips, trees[[1]])
  trees <- lapply(trees, Preorder)
  nTip <- NTip(trees[[1]])
  nTree <- length(trees)

  tr <- trees
  candidates <- character(nTip - 2L)
  score <- double(nTip - 2)
  score[1] <- ConsensusInfo(trees, info = info, check.tips = FALSE)
  cli_progress_bar("Dropping leaves", total = nTip - 2L)
  for (i in 1 + seq_len(nTip - 3L)) {
    cli_progress_update(1, status = "Remove tip {i}/{nTip - 2L}")
    tipScores <- TipInstability(tr)
    candidate <- which.max(tipScores)
    if (length(candidate)) {
      candidates[i] <- names(candidate)
    }
    tr <- lapply(tr, DropTip, candidate)
    score[i] <- ConsensusInfo(tr, info = info, check.tips = FALSE)
  }
  cli_progress_done()
  dropped <- if (fullSeq) {
    candidates
  } else {
    candidates[seq_len(which.max(score))[-1]]
  }
  score <- score[seq_len(length(dropped) + 1L)]
  # Return:
  data.frame(num = c(seq_along(score) - 1L),
             taxNum = c(NA, match(dropped, trees[[1]]$tip.label)),
             taxon = c(NA, dropped),
             rawImprovement = c(NA, score[-1] - score[-length(score)]),
             IC = score)
}

#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' cli_alert_success
#' @importFrom TreeDist ConsensusInfo
#' @importFrom TreeTools DropTip SplitFrequency Preorder RenumberTips
#' @importFrom utils combn
Roguehalla <- function (trees, dropsetSize = 1, info = 'phylogenetic') {
  if (!inherits(trees, 'multiPhylo')) {
    if (inherits(trees, 'phylo')) {
      return (trees)
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
  best <- ConsensusInfo(trees, info = info, check.tips = FALSE)

  .Drop <- function (n) {
    cli_progress_bar(paste0("Dropset size ", n))
    drops <- combn(NTip(trees[[1]]), n)
    cli_progress_update(set = 0, total = ncol(drops))
    candidates <- apply(drops, 2, function (drop) {
      cli_progress_update(1, .envir = parent.frame(2), status = paste0(
        "Drop ", startTip - NTip(trees[[1]]), " leaves = ",
        signif(best), " bits."))
      dropForest <- lapply(trees, DropTip, drop)
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
        taxSeq <- c(taxSeq, paste0(match(thisDrop, labels), collapse = ','))
        dropInf <- c(dropInf, best)
        trees <- lapply(trees, DropTip, thisDrop)
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
  data.frame(num = c(NA, seq_along(dropSeq) - 1L),
             taxNum = c(NA, taxSeq),
             taxon = c(NA, dropSeq),
             rawImprovement = c(NA, dropInf[-1] - dropInf[-length(dropInf)]),
             IC = dropInf)
}
