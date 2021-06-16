Delete <- function (f) if (file.exists(f)) file.remove(f)

test_that("Wrapper doesn't explode RogueNaRok()", {
  
  bootTrees <- system.file('example/150.bs', package = 'RogueTaxa')
  
  Rogues(ape::read.tree(bootTrees))
  skip_if(expect_true(T))
  
  
  trees <- ape::read.tree(bootTrees)
  #treeFile <- system.file('example/150.tr', package = 'RogueTaxa')
  
  neverDrop = character(0)
  computeSupport = F
  dropsetSize = 1
  labelPenalty = 0
  mreOptimization = F
  threshold = 50
  bestTree = NULL
  
  ######################
  
  wd <- tempdir()
  if (!inherits(trees, 'multiPhylo')) {
    if (inherits(trees, 'phylo')) return (NA)
    trees <- structure(trees, class = 'multiPhylo')
  }
  bootTrees <- tempfile(tmpdir = wd)
  write.tree(trees, file = bootTrees)
  on.exit(unlink(bootTrees))
  
  
  if (inherits(bestTree, 'phylo')) {
    treeFile <- tempfile(tmpdir = wd)
    write.tree(bestTree, treeFile)
    on.exit(file.remove(treeFile))

  } else {
    treeFile <- ""
  }
  if (length(neverDrop)) {
    excludeFile <- tempfile(tmpdir = wd)
    write(neverDrop, excludeFile)
    on.exit(file.remove(excludeFile))
  } else {
    excludeFile <- ""
  }
  C_RogueNaRok(bootTrees = bootTrees, 
               runId = "tmp",
               treeFile = treeFile,
               computeSupport = computeSupport,
               dropsetSize = dropsetSize,
               excludeFile = excludeFile,
               workDir = wd,
               labelPenalty = labelPenalty,
               mreOptimization = mreOptimization,
               threshold = threshold)
  rogueFile <- paste0(wd, '/RogueNaRokR_droppedRogues.tmp')
  if (!file.exists(rogueFile)) stop(rogueFile, ' not there')
  droppedRogues <- read.table(rogueFile, header = TRUE)
  unlink(rogueFile)
  unlink(paste0(wd, '/RogueNaRokR_info.tmp'))
  C_RogueNaRok(bootTrees = bootTrees, 
               runId = "tmp2",
               treeFile = treeFile,
               computeSupport = computeSupport,
               dropsetSize = dropsetSize,
               excludeFile = excludeFile,
               workDir = wd,
               labelPenalty = labelPenalty,
               mreOptimization = mreOptimization,
               threshold = threshold)
  
  
  rogueFile <- paste0(wd, '/RogueNaRokR_droppedRogues.tmp2')
  unlink(paste0(wd, '/RogueNaRokR_info.tmp2'))
  droppedRogues
  expect_true(T)
  
  
  
  Delete('RogueNaRokR_droppedRogues.tmp')
  Delete('RogueNaRokR_info.tmp')
  expect_equal(0, C_RogueNaRok(bootTrees = bootTrees,# treeFile = treeFile,
                               dropsetSize = 1,
                               runId = 'tmp'))
  
  skip_if(T)
  
  # 
  # dims <- dim(read.table('RogueNaRokR_droppedRogues.tmp', header = TRUE))
  # expect_lt(2, dims[1])
  # expect_equal(5, dims[2])
  # # if run_id exists, won't run.
  # Delete('RogueNaRokR_droppedRogues.tmp')
  # Delete('RogueNaRokR_info.tmp')
  
  Rogues(ape::read.tree(bootTrees))
  
})

test_that("C_RogueNaRok() doesn't explode", {
  skip_if(T)
  bootTrees <- system.file('example/150.bs', package = 'RogueTaxa')
  treeFile <- system.file('example/150.tr', package = 'RogueTaxa')
  # Delete('RogueNaRokR_droppedRogues.tmp')
  # Delete('RogueNaRokR_info.tmp')
  # expect_equal(0, C_RogueNaRok(bootTrees = bootTrees,# treeFile = treeFile,
  #                              dropsetSize = 1,
  #                              runId = 'tmp'))
  # 
  # dims <- dim(read.table('RogueNaRokR_droppedRogues.tmp', header = TRUE))
  # expect_lt(2, dims[1])
  # expect_equal(5, dims[2])
  # # if run_id exists, won't run.
  # Delete('RogueNaRokR_droppedRogues.tmp')
  # Delete('RogueNaRokR_info.tmp')
  
  Rogues(ape::read.tree(bootTrees))
  
})

test_that("Rogues found", {
  skip_if(T || 'Rogues found')
  library("TreeTools", warn.conflicts = FALSE, quietly = TRUE)
  trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')
  #trees <- lapply(trees, unroot)
  if (!inherits(trees, 'multiPhylo')) {
    if (inherits(trees, 'phylo')) return (NA)
    trees <- structure(trees, class = 'multiPhylo')
  }
  
  bootFile <- tempfile()
  write.tree(structure(trees, class = 'multiPhylo'), bootFile)
  on.exit(unlink(bootFile))
  Delete('RogueNaRokR_droppedRogues.tmp')
  Delete('RogueNaRokR_info.tmp')
  expect_equal(0, C_RogueNaRok(bootTrees = bootFile,# treeFile = treeFile,
                               dropsetSize = 1,
                               runId = 'tmp'))
  Delete('RogueNaRokR_droppedRogues.tmp')
  Delete('RogueNaRokR_info.tmp')
  
  expect_equal('Rogue', Rogues(trees[2:13], dropsetSize = 1L)[2, 'taxon'])
  
  skip_if(T)
  
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
