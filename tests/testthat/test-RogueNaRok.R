test_that("RogueNaRok() doesn't explode", {
  bootTrees <- system.file('example/150.bs', package = 'RogueTaxa')
  treeFile <- system.file('example/150.tr', package = 'RogueTaxa')
  expect_equal(0, RogueNaRok(bootTrees = bootTrees, treeFile = treeFile,
                             run_id = 1))
  
test_that("Rogues found", {
  library("TreeTools", warn.conflicts = FALSE, quietly = TRUE)
  trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')
  
  #expect_equal('Rogue', RogueNaRok(trees[1:3])[2, 'taxon'])
  
  #expect_equal('Rogue', RogueNaRok(trees)[2, 'taxon'])
  
  RogueNaRok(trees)
  trees[] <- lapply(trees, AddTip, 'Rogue', 'Rogue2')
  
  # Interesting aside: Majority rule consensus favours balanced splits!
  bc <- RogueNaRok(trees)
  expect_equal(1, nrow(bc))
  
  bc <- RogueNaRok(trees[-11])
  expect_equal(3, nrow(bc))
})
})
