test_that("C_RogueNaRok() runs example files", {
  Delete <- function(f) {
    f <- paste0(tmpDir,"/", f)
    if (file.exists(f)) file.remove(f)
  }

  set.seed(0)
  tmpDir <- tempdir()
  bootTrees <- system.file("example/150.bs", package = "Rogue")
  capture.output(cOutput <- C_RogueNaRok(bootTrees = bootTrees,
                                         dropsetSize = 1,
                                         labelPenalty = 0,
                                         workDir = tmpDir,
                                         runId = "tmp"))
  expect_equal(0, cOutput)

  dims <- dim(read.table(paste0(tmpDir, "/RogueNaRok_droppedRogues.tmp"),
                                header = TRUE))
  expect_lt(2, dims[1])
  expect_equal(5, dims[2])

  expect_true(Delete("RogueNaRok_droppedRogues.tmp"))
  expect_true(Delete("RogueNaRok_info.tmp"))
  set.seed(0)
  # Use just first 200 trees for faster results; will mean fewer dims though.
  expect_equal(dim(RogueTaxa(ape::read.tree(bootTrees)[1:200],
                             info = "rbic", return = "tips",
                             labelPenalty = 0, verbose = FALSE)),
               c(22, 5), tolerance = 2/28)


  treeFile <- system.file("example/150.tr", package = "Rogue")
  capture.output(cOutput <- C_RogueNaRok(bootTrees = bootTrees,
                                         treeFile = treeFile,
                                         dropsetSize = 1,
                                         labelPenalty = 0,
                                         workDir = tmpDir,
                                         runId = "tmp"))

  dims <- dim(read.table(paste0(tmpDir, "/RogueNaRok_droppedRogues.tmp"),
                         header = TRUE))
  expect_lt(2, dims[1])
  expect_equal(5, dims[2])

  Delete("RogueNaRok_droppedRogues.tmp")
  Delete("RogueNaRok_info.tmp")
})
