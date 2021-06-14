test_that("RogueNaRok() doesn't explode", {
  bootTrees <- system.file('example/150.bs', package = 'RogueTaxa')
  treeFile <- system.file('example/150.tr', package = 'RogueTaxa')
  expect_equal(0, RogueNaRok(bootTrees = bootTrees, treeFile = treeFile,
                             run_id = 1))
  
})
