test_that("RogueNaRok() doesn't explode", {
  treeFile <- system.file('example/150.bs', package = 'RogueTaxa')
  RogueNaRok(treeFile = treeFile)
})
