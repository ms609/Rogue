#' Cophenetic distance between leaves of unweighted tree
#'
#' @param x Object of class `phylo`.
#' @param nTip Integer specifying number of leaves.
#' @return `Cophenetic()` returns an unnamed integer matrix describing the
#' number of edges between each pair of edges.
#' @author Martin R. Smith, modifying algorithm by Emmanuel Paradis
#' in `ape::dist.nodes()`.
#' @importFrom ape dist.nodes
#' @keywords internal
#' @examples
#' Cophenetic(TreeTools::BalancedTree(5))
#' @useDynLib Rogue, .registration = TRUE
#' @export
Cophenetic <- function (x, nTip = length(x$tip.label), log = FALSE) {
  x <- Preorder(x)
  edge <- x$edge - 1L
  nNode <- x$Nnode
  ret <- matrix(.Call(if (log) "COPHENETIC_LOG" else "COPHENETIC",
                      n_tip = as.integer(nTip),
                      n_node = as.integer(nNode),
                      parent = as.integer(edge[, 1]),
                      child = as.integer(edge[, 2]),
                      n_edge = as.integer(dim(edge)[1])),
                nrow = nTip + if(log) 0 else nNode)

  # Return:
  if (log) {
    ret # Already cut to size
  } else {
    ret[seq_len(nTip), seq_len(nTip)]
  }
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
#' @param log Logical specifying whether to log-transform distances before
#' evaluation.
#' @param average Character specifying whether to use `'mean'` or `'median'`
#' tip distances to calculate leaf stability.
#' @param deviation Character specifying whether to use `'sd'` or `'mad'` to
#' calculate leaf stability.
#' @examples
#' library("TreeTools", quietly = TRUE)
#' trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')[3:6]
#' plot(consensus(trees), tip.col = ColByStability(trees))
#' instab <- TipInstability(trees, log = FALSE, ave = 'mean', dev = 'mad')
#' plot(ConsensusWithout(trees, names(instab[instab > 0.2])))
#' @template MRS
#' @family tip instability functions
#' @importFrom matrixStats rowMedians
#' @importFrom Rfast rowmeans rowMads rowVars
#' @export
TipInstability <- function (trees, log = TRUE, average = 'mean',
                            deviation = 'sd') {
  dists <- .TipDistances(trees, log = log)
  dims <- dim(dists)
  nTip <- dims[1]
  nTree <- dims[3]

  dists <- matrix(dists, nTip * nTip, nTree)

  whichDev <- pmatch(tolower(deviation), c('sd', 'mad'))
  if (is.na(whichDev)) {
    stop("`deviation` must be 'sd' or 'mad'")
  }
  devs <- matrix(switch(whichDev,
                        rowVars(dists, std = TRUE, parallel = TRUE),
                        rowMads(dists, parallel = TRUE)),
                 nTip, nTip)
  diag(devs) <- 0 # Faster than setting to NA, then using rowMeans(rm.na = TRUE)

  whichAve <- pmatch(tolower(average), c('mean', 'median'))
  if (is.na(whichAve)) {
    stop("`deviation` must be 'sd' or 'mad'")
  }
  aves <- matrix(switch(whichAve, rowmeans, rowMedians)(dists), nTip, nTip)

  relDevs <- devs / mean(aves[lower.tri(aves)])

  setNames(Rfast::rowmeans(relDevs), TipLabels(trees[[1]]))
}

.TipDistances <- function (trees, log = FALSE) {
  nTip <- NTip(trees)
  if (length(unique(nTip)) > 1) {
    stop("Trees must have same number of leaves")
  }
  nTip <- nTip[1]
  labels <- trees[[1]]$tip.label
  trees[-1] <- lapply(trees[-1], RenumberTips, labels)
  dists <- vapply(trees, Cophenetic, matrix(0, nTip, nTip),
                  nTip = nTip, log = log)

}

#' `ColByStability()` returns a colour reflecting the instability of each leaf.
#' @rdname TipInstability
#' @importFrom Rfast Log
#' @importFrom grDevices hcl
#' @importFrom stats cmdscale setNames
#' @importFrom TreeTools TipLabels
#' @export
ColByStability <- function (trees, log = TRUE,
                            average = 'mean', deviation = 'sd') {
  score <- TipInstability(trees, log = log, average = average,
                          deviation = deviation)
  score <- score - min(score)
  score <- score / max(score)

  # Return:
  setNames(hcl.colors(131, 'inferno')[1 + (score * 100)],
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
