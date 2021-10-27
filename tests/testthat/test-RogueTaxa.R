library("TreeTools", warn.conflicts = FALSE, quietly = TRUE)

test_that("RogueTaxa() handles bad input", {
  sameNamed <- BalancedTree(c(letters[c(1:5, 5, 6:7)]))
  bal8 <- BalancedTree(8)
  expect_error(RogueTaxa(c(sameNamed, sameNamed)))

  expect_error(RogueTaxa(c(bal8, BalancedTree(9), PectinateTree(8))))
  expect_error(RogueTaxa(c(bal8, PectinateTree(1:8))))

  expect_equal(RogueTaxa(c(bal8, bal8)), RogueTaxa(bal8))
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


  trees <- lapply(trees, AddTip, 'Rogue', 'Rogue2')

  # Interesting aside: Majority rule consensus favours balanced splits!
  bc <- RogueTaxa(trees, info = 'rbic',
                  labelPenalty = 0, verbose = FALSE)
  expect_equal(1, nrow(bc))

  bc <- RogueTaxa(trees[-11], info = 'rb',
                  labelPenalty = 0, verbose = FALSE, dropset = 2)
  expect_equal(2, nrow(bc)) # Row 1 contains a 2-taxon dropset.

  trees <- read.tree(system.file('example/150.bs', package = 'Rogue'))[1:50]
  expect_lt(1, nrow(RogueTaxa(trees, info = 'R', mreOptimization = TRUE)))
  expect_lt(1, nrow(RogueTaxa(trees, info = 'rBi', threshold = 100)))
})

test_that("Wilkinson & Crotti's examples are satisfied", {
  scaffold <- BalancedTree(c(6:4, 1:3))
  fig2 <- list(AddTip(scaffold, '3', 'X'),
               AddTip(scaffold, '4', 'X'))
  trees <- fig2
  expect_equal("X", RogueTaxa(fig2, info = 'rbic', verbose = FALSE)[2, 'taxon'])
  expect_equal("X", RogueTaxa(fig2)[2, 'taxon'])

  fig2b <- fig2[rep(1:2, c(67, 33))]
  expect_equal(1, nrow(RogueTaxa(fig2b, info = 'rbic',
                                 labelPenalty = 0, verbose = FALSE)))
  expect_equal(1, nrow(RogueTaxa(fig2b)))

  fig3 <- lapply(list(AddTip(scaffold, '1', 'X'),
                      AddTip(scaffold, '6', 'X')), AddTip, 'X', 'Y')

  trees <- fig3
  expect_equal(1, nrow(RogueTaxa(fig3, info = 'rbic', verbose = FALSE)))
  expect_equal(1, nrow(RogueTaxa(fig3)))

  fig3b <- fig3[rep(1:2, c(60, 40))]
  expect_equal(1, nrow(RogueTaxa(fig3b, info = 'rbic', verbose = FALSE)))
  expect_equal(1, nrow(RogueTaxa(fig3b)))
})
