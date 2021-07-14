#' @importFrom stats cophenetic
Cophenetic <- function (x) {
        if (is.null(x$edge.length)) {
                x$edge.length <- rep_len(1, dim(x$edge)[1])
        }
        ret <- cophenetic(x)
        ret[ret < sqrt(.Machine$double.eps)] <- 0
        ret
}

#' Tip instability
#'
#' `TipInstability()` calculates the instability of each leaf in a tree.
#' Unstable leaves are likely to display roguish behaviour.
#'
#' I define the *instability* of a pair of leaves as the median absolute
#' divergence in the cophenetic distance
#' (the number of edges in the shortest path between the leaves) across all
#' trees, normalized against the mean cophenetic distance.
#' The instability of a single leaf is the mean instability of all pairs that
#' include that leaf; higher values characterise leaves whose position is more
#' variable between trees.
#'
#' @inheritParams RogueTaxa
#' @examples
#' library("TreeTools", quietly = TRUE)
#' trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')[3:6]
#' plot(consensus(trees), tip.col = ColByStability(trees))
#' instab <- TipInstability(trees)
#' plot(ConsensusWithout(trees, names(instab[instab > 0.2])))
#' @template MRS
#' @family tip instability functions
#' @importFrom matrixStats rowMedians
#' @export
TipInstability <- function (trees) {
  dists <- .TipDistances(trees)
  dims <- dim(dists)
  nTip <- dims[1]
  nTree <- dims[3]

  n <- dim(dists)[3]
  if (n == 0) stop("No trees.")

  flatDists <- matrix(dists, nTip * nTip, nTree)
  centre <- rowMedians(flatDists)
  mads <- matrix(1.4826 * rowMedians(abs(flatDists - centre)), nTip, nTip)
  diag(mads) <- NA

  means <- rowMeans(dists, dims = 2)
  relDevs <- mads / mean(means[lower.tri(means)])
  dimnames(relDevs) <- dimnames(means)

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

#' `ColByStability()` returns a colour reflecting the instability of each leaf.
#' @rdname TipInstability
#' @param luminence Numeric luminance value to pass to [hcl()].
#' @importFrom grDevices hcl
#' @importFrom stats cmdscale setNames
#' @importFrom TreeTools TipLabels
#' @export
ColByStability <- function (trees, luminence = 50) {
  dists <- .TipDistances(trees)

  means <- rowMeans(dists, dims = 2)
  devs <- apply(dists, 1:2, function(x) mad(x))
  diag(devs) <- NA
  relDevs <- devs / mean(means[lower.tri(means)])

  pc <- cmdscale(means, k = 1)
  pc <- pc - min(pc)
  pc <- pc * 340 / max(pc)

  # Return:
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
#' sb <- TipVolatility(trees)
#' sbNorm <- 1 + (99 * (sb - min(sb)) / (max(sb - min(sb))))
#' col <- hcl.colors(128, 'inferno')[sbNorm]
#' plot(consensus(trees), tip.color = col)
#' plot(ConsensusWithout(trees, names(sb[sb == max(sb)])))
#' @importFrom TreeDist PhylogeneticInfoDistance
#' @importFrom TreeTools CladisticInfo
#' @family tip instability functions
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
QuickRogue <- function (trees, info = 'phylogenetic', fullSeq = FALSE) {
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
  candidates <- character(nTip - 2L)
  score <- double(nTip - 2)
  score[1] <- ConsensusInfo(trees, info = info, check.tips = FALSE)
  nDrops <- nTip - 3L
  cli_progress_bar("Dropping leaves", total = nDrops * (nDrops + 1L) / 2 )
  for (i in 1 + seq_len(nDrops)) {
    cli_progress_update(nDrops - (i - 1),
                        status = paste0("Leaf ", i - 1, "/", nDrops))
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
#' @importFrom TreeDist ConsensusInfo
#' @importFrom TreeTools DropTip SplitFrequency Preorder RenumberTips
#' @importFrom utils combn
Roguehalla <- function (trees, dropsetSize = 1, info = 'phylogenetic') {
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
  data.frame(num = c(0, seq_along(dropSeq)),
             taxNum = c(NA_character_, taxSeq),
             taxon = c(NA_character_, dropSeq),
             rawImprovement = c(NA, dropInf[-1] - dropInf[-length(dropInf)]),
             IC = dropInf,
             stringsAsFactors = FALSE)
}
