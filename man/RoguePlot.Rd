\name{RoguePlot}
\alias{RoguePlot}
\title{Plot rogue taxon positions across trees}
\description{
  \code{RoguePlot()} visualizes where a rogue taxon appears across different trees
  in a set, displaying a consensus tree with edges colored by the frequency
  with which the rogue taxon would be attached at each position.
}
\usage{
RoguePlot(trees, rogueTaxon, p = 0.5, legend = NULL, legend.inset = 0.02,
          pal = hcl.colors(131, "viridis")[1:101], ...)
}
\arguments{
  \item{trees}{A list of phylogenetic trees (class \code{multiPhylo} or \code{list}).}
  \item{rogueTaxon}{Character string specifying the name of the rogue taxon
    to visualize.}
  \item{p}{Numeric specifying the threshold for consensus tree construction
    (default: 0.5 for majority rule consensus).}
  \item{legend}{Character string specifying legend position, or \code{NULL} to
    suppress legend. Options include "topleft", "topright", "bottomleft",
    "bottomright", etc.}
  \item{legend.inset}{Numeric specifying legend inset from plot margins.}
  \item{pal}{A vector listing a sequence of colours to be used for plotting.
    The earliest entries will be assigned to positions where the rogue appears
    least frequently.}
  \item{...}{Additional arguments passed to plotting functions.}
}
\value{
  A list with components:
  \describe{
    \item{\code{legendLabels}}{Character vector of labels for the legend}
    \item{\code{frequencies}}{Named vector of attachment frequencies for each edge}
    \item{\code{consensus}}{The consensus tree used for plotting}
  }
}
\details{
  This function creates a consensus tree from the input trees after removing
  the specified rogue taxon. It then analyzes where the rogue taxon appears
  in each of the original trees and colors the edges of the consensus tree
  based on how frequently the rogue taxon would attach at each position.
  
  Edges with higher attachment frequencies are colored with brighter/warmer
  colors from the specified palette, while edges with lower frequencies
  receive darker/cooler colors.
  
  The function is designed to work with the expected workflow from the Rogue
  package, where users identify rogue taxa and then want to visualize their
  impact on tree topology.
}
\examples{
\dontrun{
library("ape")

# Create trees with a rogue taxon in different positions
trees <- list(
  read.tree(text = "(a, (b, (c, (rogue, (d, (e, f))))));"),
  read.tree(text = "(rogue, (a, (b, (c, (d, (e, f))))));"),
  read.tree(text = "((rogue, a), (b, (c, (d, (e, f)))));")
)

# Plot rogue positions
plotted <- RoguePlot(trees, "rogue", legend = "topleft", legend.inset = 0.02)

# Add spectrum legend (requires PlotTools package)
if (require("PlotTools", quietly = TRUE)) {
  PlotTools::SpectrumLegend(
    "bottomleft",
    palette = colorRampPalette(c(par("fg"), "#009E73"), space = "Lab")(100),
    legend = plotted$legendLabels,
    cex = 0.4
  )
}
}
}
\seealso{
  \code{\link{RogueTaxa}}, \code{\link{QuickRogue}}, \code{\link{ColByStability}}
}
\keyword{phylogeny}
\keyword{visualization}