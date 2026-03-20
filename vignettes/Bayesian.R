## ----load-packages, message = FALSE-------------------------------------------
library("TreeTools") # Read and plot trees
library("Rogue") # Find rogue taxa

## ----test-connection, echo = FALSE--------------------------------------------
online <- tryCatch({
    read.csv("https://raw.githubusercontent.com/ms609/hyoliths/master/MrBayes/hyo.nex.pstat")
    TRUE
  },
  warning = function(w) invokeRestart(""),
  error = function(e) FALSE
)

## ----load-trees---------------------------------------------------------------
if (online) {
  dataFolder <- "https://raw.githubusercontent.com/ms609/hyoliths/master/MrBayes/"
  run1.t <- paste0(dataFolder, "hyo.nex.run1.t")
  # Reading 10k trees takes a second or two...
  run1Trees <- ape::read.nexus(run1.t)
  if (packageVersion('ape') <= "5.6.1") {
    # Workaround for a bug in ape, hopefully fixed in v5.6.2
    run1Trees <- structure(lapply(run1Trees, function(tr) {
      tr$tip.label <- attr(run1Trees, 'TipLabel')
      tr
    }), class = 'multiPhylo')
  }
} else {
  # If no internet connection, we can generate some example trees instead
  run1Trees <- structure(unlist(lapply(0:21, function(backbone) {
    AddTipEverywhere(ape::as.phylo(0, nTip = 12), 'ROGUE')
  }), recursive = FALSE), class = 'multiPhylo')
}

## ----discard-burnin-----------------------------------------------------------
burninFrac <- 0.25
nTrees <- length(run1Trees)
trees <- run1Trees[seq(from = burninFrac * nTrees, to = nTrees)]

## ----thin-trees---------------------------------------------------------------
sampleSize <- 100
trees <- run1Trees[seq(from = burninFrac * nTrees, to = nTrees,
                       length.out = sampleSize)]

## ----initial-view, fig.width = 8, fig.asp = 2/3-------------------------------
plenary <- Consensus(trees, p = 0.5)

par(mar = rep(0, 4), cex = 0.85)
plot(plenary, tip.color = ColByStability(trees))
PlotTools::SpectrumLegend(
  "bottomright", legend = c("Stable", "Unstable"),
  palette = hcl.colors(131, 'inferno')[1:101]
)

## ----find-rogues--------------------------------------------------------------
rogues <- QuickRogue(trees)
# rogues <- RogueTaxa(trees) might do a better job, much more slowly
rogues

# The first line reports the information content of the plenary tree
rogueTaxa <- rogues$taxon[-1]

## ----trees, fig.asp = 1.2/2, fig.width = 8------------------------------------
par(mar = rep(0, 4)) # Remove plot margin
par(mfrow = c(1, 2)) # Multi-panel plot
par(cex = 0.85) # Smaller labels

plenary <- Consensus(trees, p = 0.5)
reduced <- ConsensusWithout(trees, rogueTaxa, p = 0.5)

plot(plenary,
     tip.color = ifelse(plenary$tip.label %in% rogueTaxa, 2, 1))
LabelSplits(plenary, SplitFrequency(plenary, trees))
plot(reduced)
LabelSplits(reduced, SplitFrequency(reduced, trees))

