Delete <- function (f) if (file.exists(f)) file.remove(f)
test_that("RogueNaRok() doesn't explode", {
  bootTrees <- system.file('example/150.bs', package = 'RogueTaxa')
  treeFile <- system.file('example/150.tr', package = 'RogueTaxa')
  Delete('RogueNaRokR_droppedRogues.tmp')
  Delete('RogueNaRokR_info.tmp')
  expect_equal(0, C_RogueNaRok(bootTrees = bootTrees,# treeFile = treeFile,
                               dropsetSize = 1,
                               runId = 'tmp'))
  
  dims <- dim(read.table('RogueNaRokR_droppedRogues.tmp', header = TRUE))
  expect_lt(2, dims[1])
  expect_equal(5, dims[2])
  # if run_id exists, won't run.
  Delete('RogueNaRokR_droppedRogues.tmp')
  Delete('RogueNaRokR_info.tmp')
  
})

test_that("Rogues found", {
  library("TreeTools", warn.conflicts = FALSE, quietly = TRUE)
  trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')
  expect_equal('Rogue', Rogues(trees[2:3], dropsetSize = 1L)[2, 'taxon'])
  
  
  expect_equal('Rogue', Rogues(trees)[2, 'taxon'])
  
  Rogues(trees)
  trees[] <- lapply(trees, AddTip, 'Rogue', 'Rogue2')
  
  # Interesting aside: Majority rule consensus favours balanced splits!
  bc <- Rogues(trees)
  expect_equal(1, nrow(bc))
  
  bc <- Rogues(trees[-11])
  expect_equal(3, nrow(bc))
})

test_that("Wilkinson & Crotti's examples are satisfied", {
  skip_if(TRUE)
  scaffold <- BalancedTree(c(6:4, 1:3))
  fig2 <- list(AddTip(scaffold, '3', 'X'),
               AddTip(scaffold, '4', 'X'))
  trees <- fig2
  expect_equal(Rogues(fig2)[2, 'taxon'])
  
  fig2b <- fig2[rep(1:2, c(67, 33))]
  expect_equal(Rogues(fig2b)[2, 'taxon'])
  
  fig3 <- lapply(list(AddTip(scaffold, '1', 'X'),
                      AddTip(scaffold, '6', 'X')), AddTip, 'X', 'Y')
  
  trees <- fig3
  expect_equal(c('X', 'Y'), Rogues(fig3)[2:3, 'taxon'])
  
  fig3b <- fig3[rep(1:2, c(60, 40))]
  expect_equal(c('X', 'Y'), Rogues(fig3b)[2:3, 'taxon'])
  
  fig3c <- lapply(fig3b, drop.tip, names(tr3b[tr3b == max(tr3b)])) 
  expect_true(all(TipVolatility(fig3c) == 0))
  expect_equal(1, nrow(Rogues(fig3b)))
})
