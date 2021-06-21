Delete <- function (f) if (file.exists(f)) file.remove(f)

test_that("C_RogueNaRok() doesn't explode", {
  set.seed(0)
  bootTrees <- system.file('example/150.bs', package = 'RogueTaxa')
  treeFile <- system.file('example/150.tr', package = 'RogueTaxa')
  Delete('RogueNaRokR_droppedRogues.tmp')
  Delete('RogueNaRokR_info.tmp')
  capture.output(cOutput <- C_RogueNaRok(bootTrees = bootTrees,# treeFile = treeFile,
                                         dropsetSize = 1,
                                         labelPenalty = 0,
                                         runId = 'tmp'))
  expect_equal(0, cOutput)

  dims <- dim(read.table('RogueNaRokR_droppedRogues.tmp', header = TRUE))
  expect_lt(2, dims[1])
  expect_equal(5, dims[2])
  # if run_id exists, won't run.
  Delete('RogueNaRokR_droppedRogues.tmp')
  Delete('RogueNaRokR_info.tmp')
  
  expect_equal(dims, dim(RogueTaxa(ape::read.tree(bootTrees), 
                                   labelPenalty = 0, verbose = FALSE)),
               tolerance = 2/28)
})

test_that("Rogues found", {
  library("TreeTools", warn.conflicts = FALSE, quietly = TRUE)
  trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')
  if (!inherits(trees, 'multiPhylo')) {
    if (inherits(trees, 'phylo')) return (NA)
    trees <- structure(trees, class = 'multiPhylo')
  }
  
  expect_equal('Rogue', RogueTaxa(trees[2:13], dropsetSize = 1L,
                                  labelPenalty = 0,
                                  verbose = FALSE)[2, 'taxon'])
  expect_equal('Rogue', RogueTaxa(trees, labelPenalty = 0,
                                  verbose = FALSE)[2, 'taxon'])
  
  
  trees[] <- lapply(trees, AddTip, 'Rogue', 'Rogue2')
  
  # Interesting aside: Majority rule consensus favours balanced splits!
  bc <- RogueTaxa(trees,
                  labelPenalty = 0, verbose = FALSE)
  expect_equal(1, nrow(bc))
  
  bc <- RogueTaxa(trees[-11],
                  labelPenalty = 0, verbose = FALSE, dropset = 2)
  expect_equal(2, nrow(bc)) # Row 1 contains a 2-taxon dropset.
})

test_that("Wilkinson & Crotti's examples are satisfied", {
  scaffold <- BalancedTree(c(6:4, 1:3))
  fig2 <- list(AddTip(scaffold, '3', 'X'),
               AddTip(scaffold, '4', 'X'))
  trees <- fig2
  expect_equal("X", RogueTaxa(fig2, verbose = FALSE)[2, 'taxon'])
  
  skip_if(TRUE)
  fig2b <- fig2[rep(1:2, c(67, 33))]
  expect_equal('X', RogueTaxa(fig2b, verbose = FALSE)[2, 'taxon'])
  expect_equal(NA, RogueTaxa(fig2b, labelPenalty = 0, verbose = FALSE)[2, 'taxon'])
  
  fig3 <- lapply(list(AddTip(scaffold, '1', 'X'),
                      AddTip(scaffold, '6', 'X')), AddTip, 'X', 'Y')
  
  trees <- fig3
  expect_equal(c('X', 'Y'), RogueTaxa(fig3, verbose = FALSE)[2:3, 'taxon'])
  
  fig3b <- fig3[rep(1:2, c(60, 40))]
  expect_equal(c('X', 'Y'), RogueTaxa(fig3b, verbose = FALSE)[2:3, 'taxon'])
  
  fig3c <- lapply(fig3b, drop.tip, names(tr3b[tr3b == max(tr3b)])) 
  expect_true(all(TipVolatility(fig3c) == 0))
  expect_equal(1, nrow(RogueTaxa(fig3b)))
})
