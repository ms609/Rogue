#' Plot rogue taxon positions across trees
#'
#' `RoguePlot()` visualizes where a rogue taxon appears across different trees
#' in a set, displaying a consensus tree with edges colored by the frequency
#' with which the rogue taxon would be attached at each position.
#'
#' @param trees A list of phylogenetic trees (class `multiPhylo` or `list`).
#' @param rogueTaxon Character string specifying the name of the rogue taxon
#'   to visualize.
#' @param p Numeric specifying the threshold for consensus tree construction
#'   (default: 0.5 for majority rule consensus).
#' @param legend Character string specifying legend position, or `NULL` to
#'   suppress legend. Options include "topleft", "topright", "bottomleft",
#'   "bottomright", etc.
#' @param legend.inset Numeric specifying legend inset from plot margins.
#' @param pal A vector listing a sequence of colours to be used for plotting.
#'   The earliest entries will be assigned to positions where the rogue appears
#'   least frequently.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return A list with components:
#'   \describe{
#'     \item{`legendLabels`}{Character vector of labels for the legend}
#'     \item{`frequencies`}{Named vector of attachment frequencies for each edge}
#'     \item{`consensus`}{The consensus tree used for plotting}
#'   }
#'
#' @examples
#' \dontrun{
#' library("TreeTools", quietly = TRUE)
#'
#' # Create trees with a rogue taxon in different positions
#' trees <- list(
#'   read.tree(text = "(a, (b, (c, (rogue, (d, (e, f))))));"),
#'   read.tree(text = "(rogue, (a, (b, (c, (d, (e, f))))));"),
#'   read.tree(text = "((rogue, a), (b, (c, (d, (e, f)))));")
#' )
#' 
#' # Plot rogue positions
#' plotted <- RoguePlot(trees, "rogue", legend = "topleft", legend.inset = 0.02)
#' 
#' # Add spectrum legend
#' PlotTools::SpectrumLegend(
#'   "bottomleft",
#'   palette = colorRampPalette(c(par("fg"), "#009E73"), space = "Lab")(100),
#'   legend = plotted$legendLabels,
#'   cex = 0.4
#' )
#' }
#'
#' @importFrom ape consensus drop.tip plot.phylo
#' @importFrom grDevices colorRampPalette hcl.colors
#' @importFrom graphics legend
#' @export
RoguePlot <- function(trees, rogueTaxon, p = 0.5, 
                      legend = NULL, legend.inset = 0.02,
                      pal = hcl.colors(131, "viridis")[1:101], ...) {
  
  # Validate input
  if (missing(trees) || missing(rogueTaxon)) {
    stop("Both 'trees' and 'rogueTaxon' arguments are required")
  }
  
  if (!inherits(trees, c("list", "multiPhylo"))) {
    stop("'trees' must be a list or multiPhylo object")
  }
  
  if (length(trees) < 2) {
    stop("At least two trees are required")
  }
  
  # Check if rogue taxon exists in trees
  allTips <- unique(unlist(lapply(trees, function(tree) tree$tip.label)))
  if (!rogueTaxon %in% allTips) {
    stop("Rogue taxon '", rogueTaxon, "' not found in trees")
  }
  
  # Remove trees that don't contain the rogue taxon
  hasRogue <- sapply(trees, function(tree) rogueTaxon %in% tree$tip.label)
  if (!any(hasRogue)) {
    stop("Rogue taxon not found in any trees")
  }
  trees <- trees[hasRogue]
  
  # Create consensus without rogue taxon
  treesWithoutRogue <- lapply(trees, function(tree) {
    if (rogueTaxon %in% tree$tip.label) {
      drop.tip(tree, rogueTaxon)
    } else {
      tree
    }
  })
  
  # Convert to multiPhylo class if needed
  class(treesWithoutRogue) <- "multiPhylo"
  
  consensus <- consensus(treesWithoutRogue, p = p)
  nEdges <- nrow(consensus$edge)
  
  # Calculate attachment frequencies for each possible position
  # This algorithm examines each tree to determine where the rogue would
  # most likely attach to the consensus topology
  frequencies <- rep(0, nEdges)
  names(frequencies) <- paste0("edge_", seq_len(nEdges))
  
  for (i in seq_along(trees)) {
    tree <- trees[[i]]
    if (!rogueTaxon %in% tree$tip.label) next
    
    # Find sister taxa of the rogue in this tree
    rogueTipIndex <- which(tree$tip.label == rogueTaxon)
    if (length(rogueTipIndex) == 0) next
    
    # Get the edge connected to the rogue tip
    rogueEdgeIndex <- which(tree$edge[, 2] == rogueTipIndex)
    if (length(rogueEdgeIndex) == 0) next
    
    rogueParentNode <- tree$edge[rogueEdgeIndex, 1]
    
    # Find sister edges (other children of the same parent)
    sisterEdges <- which(tree$edge[, 1] == rogueParentNode & tree$edge[, 2] != rogueTipIndex)
    
    # Get taxa that are sisters to the rogue
    sisterTaxa <- c()
    for (sEdge in sisterEdges) {
      descendantNode <- tree$edge[sEdge, 2]
      if (descendantNode <= length(tree$tip.label)) {
        # It's a tip
        sisterTaxa <- c(sisterTaxa, tree$tip.label[descendantNode])
      } else {
        # It's an internal node - get all descendant tips
        descendants <- which(tree$edge[, 1] == descendantNode)
        tipQueue <- tree$edge[descendants, 2]
        while (length(tipQueue) > 0) {
          currentNode <- tipQueue[1]
          tipQueue <- tipQueue[-1]
          if (currentNode <= length(tree$tip.label)) {
            sisterTaxa <- c(sisterTaxa, tree$tip.label[currentNode])
          } else {
            newDescendants <- which(tree$edge[, 1] == currentNode)
            tipQueue <- c(tipQueue, tree$edge[newDescendants, 2])
          }
        }
      }
    }
    
    # Now find corresponding edges in consensus tree
    # This is simplified - we assign frequency to edges that connect
    # to any of the sister taxa
    for (j in seq_len(nEdges)) {
      edgeDescendant <- consensus$edge[j, 2]
      if (edgeDescendant <= length(consensus$tip.label)) {
        # Edge leads to a tip
        if (consensus$tip.label[edgeDescendant] %in% sisterTaxa) {
          frequencies[j] <- frequencies[j] + 1
        }
      } else {
        # Edge leads to internal node - check if any sister taxa are descendants
        # Get all descendant tips of this internal node
        descendants <- c()
        tipQueue <- edgeDescendant
        while (length(tipQueue) > 0) {
          currentNode <- tipQueue[1]
          tipQueue <- tipQueue[-1]
          childEdges <- which(consensus$edge[, 1] == currentNode)
          for (cEdge in childEdges) {
            childNode <- consensus$edge[cEdge, 2]
            if (childNode <= length(consensus$tip.label)) {
              descendants <- c(descendants, consensus$tip.label[childNode])
            } else {
              tipQueue <- c(tipQueue, childNode)
            }
          }
        }
        
        # Check if any sister taxa are in descendants
        if (any(sisterTaxa %in% descendants)) {
          frequencies[j] <- frequencies[j] + 1
        }
      }
    }
  }
  
  # Add base frequency to all edges to avoid completely black edges
  frequencies <- frequencies + 0.1
  
  # Normalize frequencies to [0, 1] range
  if (max(frequencies) > 0) {
    frequencies <- frequencies / max(frequencies)
  }
  
  # Map frequencies to colors
  colorIndices <- 1 + as.integer(frequencies * (length(pal) - 1))
  edgeColors <- pal[colorIndices]
  
  # Create plot
  plot(consensus, edge.color = edgeColors, ...)
  
  # Add legend if requested
  if (!is.null(legend)) {
    legendLabels <- c("Low frequency", "", "High frequency")
    legend(legend, inset = legend.inset, 
           legend = legendLabels,
           fill = pal[c(1, length(pal) %/% 2, length(pal))],
           title = paste("Attachment frequency"))
  }
  
  # Create legend labels for external use
  freqRange <- range(frequencies)
  legendLabels <- sprintf("%.1f", seq(freqRange[1], freqRange[2], length.out = 4))
  
  # Return result
  invisible(list(
    legendLabels = legendLabels,
    frequencies = frequencies,
    consensus = consensus
  ))
}