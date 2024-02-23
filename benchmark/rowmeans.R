#' library("TreeTools", quietly = TRUE)
#' 
#' # Generate some trees with a rogue taxon
#' trees <- AddTipEverywhere(BalancedTree(8), "Rogue")[3:6]
#' 
#' # Display the strict consensus
#' plot(consensus(trees), tip.col = ColByStability(trees))
#' 
#' # Add a legend for the colour scale used
#' PlotTools::SpectrumLegend(
#'   "bottomleft", bty = "n", # No box
#'   legend = c("Unstable", "", "Stable"),
#'   palette = hcl.colors(131, "inferno")[1:101]
#' )
#' 
#' # Calculate leaf stability
#' instab <- TipInstability(trees, log = FALSE, ave = "mean", dev = "mad")
#' 
#' # Plot a consensus that omits the least stable leaves
#' plot(ConsensusWithout(trees, names(instab[instab > 0.2])))
#' 
#' function(trees, 
trees <- lapply(1:550, function(X) AddTip(BalancedTree(56), label = "Rogue"))

ub(TipInstability(trees, parallel = TRUE),
   TipInstability(trees, parallel = FALSE)) # Makes basically no difference

ub(TipInstability(trees, deviation = "mad", parallel = TRUE),
   TipInstability(trees, deviation = "mad", parallel = FALSE),
   times = 25) # parallel maybe 10% faster

log = TRUE
average = "mean"
deviation = "sd"
checkTips = TRUE

if (inherits(trees, "phylo") || length(trees) < 2L) {
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
                       class = "multiPhylo"))
  nTip <- nTip[1]
} else {
  nTip <- NTip(trees[[1]])
}

dists <- vapply(trees, GraphGeodesic, double(nTip * nTip),
                nTip = nTip, log = log, asMatrix = FALSE)

whichDev <- pmatch(tolower(deviation), c("sd", "mad"))
if (is.na(whichDev)) {
  stop("`deviation` must be 'sd' or 'mad'")
}
devs <- matrix(switch(whichDev,
                      rowVars(dists, std = TRUE, parallel = TRUE),
                      rowMads(dists, parallel = TRUE)),
               nTip, nTip)
devs[is.nan(devs)] <- 0 # rowVars returns NaN instead of 0
#diag(devs) <- 0 # Faster than setting to NA, then using rowMeans(rm.na = TRUE)


whichAve <- pmatch(tolower(average), c("mean", "median"))
if (is.na(whichAve)) {
  stop("`average` must be 'mean' or 'median'")
}
ub(Rfast::rowmeans(dists), Rfast::colmeans(dists))
aves <- matrix(switch(whichAve, rowmeans, rowMedians)(dists), nTip, nTip)

isSymmetric(aves)
relDevs <- devs / mean(aves[lower.tri(aves)])

ub <- function(a, b) {
  stopifnot(all.equal(a, b))
  microbenchmark::microbenchmark(a, b)
}
ub(Rfast::rowmeans(relDevs),
   Rfast::colmeans(relDevs, parallel = TRUE),
   Rfast::colmeans(relDevs, parallel = FALSE),
   times = 1000
   )

setNames(Rfast::rowmeans(relDevs), TipLabels(trees[[1]]))
