library("TreeTools", warn.conflicts = FALSE, quietly = TRUE)

test_that("RogueTaxa() handles bad input", {
  sameNamed <- BalancedTree(c(letters[c(1:5, 5, 6:7)]))
  expect_error(RogueTaxa(c(sameNamed, sameNamed)))

  expect_error(RogueTaxa(c(BalancedTree(8), BalancedTree(9), PectinateTree(8))))
  expect_error(RogueTaxa(c(BalancedTree(8), PectinateTree(1:8))))

})

test_that("C_RogueNaRok() runs example files", {
  Delete <- function (f) {
    f <- paste0(tmpDir,'/', f)
    if (file.exists(f)) file.remove(f)
  }

  set.seed(0)
  tmpDir <- tempdir()
  bootTrees <- system.file('example/150.bs', package = 'Rogue')
  capture.output(cOutput <- C_RogueNaRok(bootTrees = bootTrees,
                                         dropsetSize = 1,
                                         labelPenalty = 0,
                                         workDir = tmpDir,
                                         runId = 'tmp'))
  expect_equal(0, cOutput)

  dims <- dim(read.table(paste0(tmpDir, '/RogueNaRok_droppedRogues.tmp'),
                                header = TRUE))
  expect_lt(2, dims[1])
  expect_equal(5, dims[2])

  expect_true(Delete('RogueNaRok_droppedRogues.tmp'))
  expect_true(Delete('RogueNaRok_info.tmp'))
  set.seed(0)
  # Use just first 200 trees for faster results; will mean fewer dims though.
  expect_equal(c(22, 5), dim(RogueTaxa(ape::read.tree(bootTrees)[1:200],
                                   info = 'rbic', return = 'tips',
                                   labelPenalty = 0, verbose = FALSE)),
               tolerance = 2/28)


  treeFile <- system.file('example/150.tr', package = 'Rogue')
  capture.output(cOutput <- C_RogueNaRok(bootTrees = bootTrees,
                                         treeFile = treeFile,
                                         dropsetSize = 1,
                                         labelPenalty = 0,
                                         workDir = tmpDir,
                                         runId = 'tmp'))

  dims <- dim(read.table(paste0(tmpDir, '/RogueNaRok_droppedRogues.tmp'),
                         header = TRUE))
  expect_lt(2, dims[1])
  expect_equal(5, dims[2])

  Delete('RogueNaRok_droppedRogues.tmp')
  Delete('RogueNaRok_info.tmp')
})

test_that("Rogues found", {
  trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')
  if (!inherits(trees, 'multiPhylo')) {
    if (inherits(trees, 'phylo')) return (NA)
    trees <- structure(trees, class = 'multiPhylo')
  }

  expect_equal('Rogue', RogueTaxa(trees[2:13], info = 'rbic', dropsetSize = 1L,
                                  labelPenalty = 0,
                                  verbose = FALSE)[2, 'taxon'])
  expect_equal('Rogue', RogueTaxa(trees, info = 'rbic', labelPenalty = 0,
                                  verbose = FALSE)[2, 'taxon'])


  trees[] <- lapply(trees, AddTip, 'Rogue', 'Rogue2')

  # Interesting aside: Majority rule consensus favours balanced splits!
  bc <- RogueTaxa(trees, info = 'rbic',
                  labelPenalty = 0, verbose = FALSE)
  expect_equal(1, nrow(bc))

  bc <- RogueTaxa(trees[-11], info = 'rbic',
                  labelPenalty = 0, verbose = FALSE, dropset = 2)
  expect_equal(2, nrow(bc)) # Row 1 contains a 2-taxon dropset.

  trees <- read.tree(system.file('example/150.bs', package = 'Rogue'))[1:50]
  expect_lt(1, nrow(RogueTaxa(trees, info = 'rbic', mreOptimization = TRUE)))
  expect_lt(1, nrow(RogueTaxa(trees, info = 'rbic', threshold = 100)))
})

test_that("Wilkinson & Crotti's examples are satisfied", {
  scaffold <- BalancedTree(c(6:4, 1:3))
  fig2 <- list(AddTip(scaffold, '3', 'X'),
               AddTip(scaffold, '4', 'X'))
  trees <- fig2
  expect_equal("X", RogueTaxa(fig2, verbose = FALSE)[2, 'taxon'])

  fig2b <- fig2[rep(1:2, c(67, 33))]
  expect_equal(NA, RogueTaxa(fig2b, labelPenalty = 0, verbose = FALSE)[2, 'taxon'])

  fig3 <- lapply(list(AddTip(scaffold, '1', 'X'),
                      AddTip(scaffold, '6', 'X')), AddTip, 'X', 'Y')

  trees <- fig3
  expect_equal(1, nrow(RogueTaxa(fig3, verbose = FALSE)))

  fig3b <- fig3[rep(1:2, c(60, 40))]
  expect_equal(1, nrow(RogueTaxa(fig3b, verbose = FALSE)))
})
