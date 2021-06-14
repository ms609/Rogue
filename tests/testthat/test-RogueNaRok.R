test_that("RogueNaRok() doesn't explode", {
  bootTrees <- system.file('example/150.bs', package = 'RogueTaxa')
  treeFile <- system.file('example/150.tr', package = 'RogueTaxa')
  RogueNaRok(bootTrees = bootTrees, treeFile = treeFile)
  
})
