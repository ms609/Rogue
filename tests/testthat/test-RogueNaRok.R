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
