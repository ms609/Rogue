#' Graph Geodesic between leaves of unweighted tree
#'
#' @param x Object of class `phylo`.
#' @param nTip Integer specifying number of leaves.
#' @param asMatrix Logical specifying whether to coerce output to matrix format.
#' @return `GraphGeodesic()` returns an unnamed integer matrix describing the
#' number of edges between each pair of edges.
#' @author Martin R. Smith, modifying algorithm by Emmanuel Paradis
#' in `ape::dist.nodes()`.
#' @importFrom ape dist.nodes
#' @keywords internal
#' @examples
#' GraphGeodesic(TreeTools::BalancedTree(5))
#' @useDynLib Rogue, .registration = TRUE
#' @export
GraphGeodesic <- function (x, nTip = length(x$tip.label), log = FALSE,
                           asMatrix = TRUE) {
  x <- Preorder(x)
  edge <- x$edge - 1L
  nNode <- x$Nnode
  ret <- .Call(if (isTRUE(log)) "LOG_GRAPH_GEODESIC" else "GRAPH_GEODESIC",
               n_tip = as.integer(nTip),
               n_node = as.integer(nNode),
               parent = as.integer(edge[, 1]),
               child = as.integer(edge[, 2]),
               n_edge = as.integer(dim(edge)[1]))

  # Return:
  if (log) {
    if (asMatrix) {
      matrix(ret, nTip)
    } else {
      ret
    }
  } else {
    matrix(ret, nrow = nTip + nNode)[seq_len(nTip), seq_len(nTip)]
  }
}

#' @rdname GraphGeodesic
#' @export
Cophenetic <- function (x, nTip = length(x$tip.label), log = FALSE,
                        asMatrix = TRUE) {
  .Deprecated('GraphGeodesic')
  GraphGeodesic(x, nTip, log, asMatrix)
}

#' Tip instability
#'
#' `TipInstability()` calculates the instability of each leaf in a tree.
#' Unstable leaves are likely to display roguish behaviour.
#'
#' \insertCite{SmithCons;textual}{Rogue} defines the *instability* of a pair
#' of leaves as the median absolute divergence in the graph geodesic
#' (the number of edges in the shortest path between the leaves) across all
#' trees, normalized against the mean graph geodesic.
#' The instability of a single leaf is the mean instability of all pairs that
#' include that leaf; higher values characterise leaves whose position is more
#' variable between trees.
#'
#' Other concepts of leaf instability include
#'
#' - The 'taxonomic instability index', as implemented in Mesquite:
#' described by \insertCite{Thomson2010;textual}{Rogue} as
#' \eqn{\sum\limits_{(x, y), j \neq i}{\frac{|D~ijx~ - D~ijy~|}{(D~ijx~ - D~ijy~)^2}}}{\sum[x, y, j != i] (D[ijx] - D[ijy] / (D[ijx] - D[ijy])^2 )},
#' where \eqn{D~ijx~}{D[ijx]} is the patristic distance (i.e. length of edges)
#' between leaves \eqn{i} and \eqn{j} in tree \eqn{x}.
#'
#' - the average stability of triplets (i.e. quartets including the root) that
#' include the leaf \insertCite{Thorley1999}{Rogue}, implemented in "Phyutility"
#' \insertCite{Smith2008}{Rogue}.
#'
#' @inheritParams RogueTaxa
#' @param log Logical specifying whether to log-transform distances when
#' calculating leaf stability.
#' @param average Character specifying whether to use `'mean'` or `'median'`
#' tip distances to calculate leaf stability.
#' @param deviation Character specifying whether to use `'sd'` or `'mad'` to
#' calculate leaf stability.
#' @param checkTips Logical specifying whether to check that tips are numbered
#' consistently.
#' @references
#' \insertAllCited{}
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
                            deviation = 'sd',
                            checkTips = TRUE) {
  if (inherits(trees, 'phylo') || length(trees) < 2L) {
    tips <- TipLabels(trees)
    return(setNames(double(length(tips)), tips))
  }
  labels <- trees[[1]]$tip.label
  if (checkTips) {
    nTip <- NTip(trees)
    if (length(unique(nTip)) > 1) {
      stop("Trees must have same number of leaves")
    }
    trees <- c(trees[[1]],
               structure(lapply(trees[-1], RenumberTips, labels),
                         class = 'multiPhylo'))
    nTip <- nTip[1]
  } else {
    nTip <- NTip(trees[[1]])
  }
  nTree <- length(trees)

  dists <- vapply(trees, GraphGeodesic, double(nTip * nTip),
                  nTip = nTip, log = log, asMatrix = FALSE)

  whichDev <- pmatch(tolower(deviation), c('sd', 'mad'))
  if (is.na(whichDev)) {
    stop("`deviation` must be 'sd' or 'mad'")
  }
  devs <- matrix(switch(whichDev,
                        rowVars(dists, std = TRUE, parallel = TRUE),
                        rowMads(dists, parallel = TRUE)),
                 nTip, nTip)
  devs[is.nan(devs)] <- 0 # rowVars returns NaN instead of 0
  #diag(devs) <- 0 # Faster than setting to NA, then using rowMeans(rm.na = TRUE)


  whichAve <- pmatch(tolower(average), c('mean', 'median'))
  if (is.na(whichAve)) {
    stop("`average` must be 'mean' or 'median'")
  }
  aves <- matrix(switch(whichAve, rowmeans, rowMedians)(dists), nTip, nTip)

  relDevs <- devs / mean(aves[lower.tri(aves)])

  setNames(Rfast::rowmeans(relDevs), TipLabels(trees[[1]]))
}

#' `ColByStability()` returns a colour reflecting the instability of each leaf.
#' @rdname TipInstability
#' @importFrom Rfast Log
#' @importFrom grDevices hcl.colors
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
#' phylogenetic information distance \insertCite{Smith2020}{Rogue} between
#' trees when that tip is removed.
#' Effective when the number of trees is small.
#'
#' @inheritParams RogueTaxa
#' @return `TipVolatility()` returns a named vector listing the volatility
#' index calculated for each leaf.
#' @references
#' \insertAllCited{}
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
#' @importFrom TreeTools CladisticInfo DropTipPhylo
#' @family tip instability functions
#' @export
TipVolatility <- function (trees) {
  tips <- trees[[1]]$tip.label
  startInfo <- mean(CladisticInfo(trees))
  info <- vapply(tips, function (drop) {
    tr <- lapply(trees, DropTipPhylo, drop, check = FALSE)
    c(meanInfo = mean(CladisticInfo(tr)),
      meanDist = mean(PhylogeneticInfoDistance(tr, normalize = TRUE)))
  }, double(2))
  mean(PhylogeneticInfoDistance(trees, normalize = TRUE)) - info['meanDist', ]
}
