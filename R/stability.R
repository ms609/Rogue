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
Cophenetic <- function (x, nTip = length(x$tip.label)) {
  x <- Preorder(x)
  edge <- x$edge - 1L
  nNode <- x$Nnode
  ret <- matrix(.Call("COPHENETIC",
                      n_tip = as.integer(nTip),
                      n_node = as.integer(nNode),
                      parent = as.integer(edge[, 1]),
                      child = as.integer(edge[, 2]),
                      n_edge = as.integer(dim(edge)[1])), nTip + nNode)

  # Return:
  ret[seq_len(nTip), seq_len(nTip)]
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
  labels <- trees[[1]]$tip.label
  trees[-1] <- lapply(trees[-1], RenumberTips, labels)
  dists <- vapply(trees, Cophenetic,
                  matrix(0, nTip, nTip, dimnames = list(labels, labels)),
                  nTip = nTip)

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
  dims <- dim(dists)
  nTip <- dims[1]
  nTree <- dims[3]

  flatDists <- matrix(dists, nTip * nTip, nTree)
  centre <- rowMedians(flatDists)
  mads <- matrix(1.4826 * rowMedians(abs(flatDists - centre)), nTip, nTip)
  diag(mads) <- NA

  means <- rowMeans(dists, dims = 2)
  relDevs <- mads / mean(means[lower.tri(means)])

  pc <- cmdscale(means, k = 1)
  pc <- pc - min(pc)
  pc <- pc * 340 / max(pc)

  # Return:
  setNames(hcl(h = pc, c = 100 * (1 - rowMeans(relDevs, na.rm = TRUE)),
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
