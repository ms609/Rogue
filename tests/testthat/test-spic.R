test_that("Roguehalla() handles odd input", {
  trees <- list(ape::read.tree(text = '(a, (b, (c, (d, (e, X)))));'),
                ape::read.tree(text = '((a, X), (b, (c, (d, e))));'))
  ic <- ConsensusInfo(lapply(trees, DropTip, 'X'), 'p')
  expect_equal(data.frame(num = c(NA, 0),
                           taxNum = c(NA_character_, '6'),
                           taxon = c(NA_character_, 'X'),
                           rawImprovement = c(NA, ic),
                           IC = c(0, ic)),
                Roguehalla(trees))

  expect_error(RogueTaxa(trees, info = 'Error'))
  expect_error(RogueTaxa(trees, return = 'Error'))
})

test_that("Rogues found", {

  library("TreeTools", quietly = TRUE)
  trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')
  instab <- TipInstability(trees)
  expect_equal('Rogue', names(which.max(instab)))

  ci <- TipVolatility(trees)
  expect_equal('Rogue', names(which.max(ci)))

  dists <- TreeDist::PhylogeneticInfoDistance(trees, normalize = TRUE)
  expect_equal(mean(dists) - 0, max(ci))

  expect_equal(2L, nrow(QuickRogue(trees)))
  expect_equal(8L, NTip(RogueTaxa(trees, return = 'TREE')))
  expect_equal(8L, NTip(RogueTaxa(trees, info = 'fsp', return = 'tr')))
  expect_equal(2L, nrow(Roguehalla(trees, 1)))
  expect_equal(8L, NTip(RogueTaxa(trees, info = 'sp', return = 'tr')))


  trees[] <- lapply(trees, AddTip, 'Rogue', 'Rogue2')
  ci <- TipVolatility(trees)
  expect_equal(c('Rogue', 'Rogue2'), names(ci[ci == max(ci)]))

  # Interesting aside: Majority rule consensus favours balanced splits!
  bc <- RogueTaxa(trees, info = 'fsp', return = 'TR')
  expect_equal(10L, NTip(bc))
  expect_equal(9L, bc$Nnode)

  bc <- RogueTaxa(trees, info = 'fsc', return = 'TR')
  expect_equal(10L, NTip(bc))
  expect_equal(9L, bc$Nnode)

  bc <- RogueTaxa(trees[-11], info = 'fsp', return = 'TR')
  expect_equal(8L, NTip(bc))
  expect_equal(7L, bc$Nnode)

  expect_equal(1L, nrow(RogueTaxa(trees[-11], drop = 1, 'r')))
  expect_equal(1L, nrow(RogueTaxa(trees[-11], drop = 1, 'phy')))
  expect_equal(1L, nrow(RogueTaxa(trees[-11], drop = 1, 'clu')))
  expect_equal(2L, nrow(RogueTaxa(trees[-11], drop = 2, 'r')))
  expect_equal(2L, nrow(RogueTaxa(trees[-11], drop = 2, 'sp')))
  expect_equal(2L, nrow(RogueTaxa(trees[-11], drop = 2, 'sci')))
  # Check trees are created
  expect_equal(8L, NTip(RogueTaxa(trees[-11], dr = 2, ret = 'Tr', 'r')))
  expect_equal(8L, NTip(RogueTaxa(trees[-11], dr = 2, ret = 'Tr', 'sp')))
  expect_equal(8L, NTip(RogueTaxa(trees[-11], dr = 2, ret = 'Tr', 'sci')))

})
#
# test_that("Benchmarking", {
#         skip_if(TRUE)
#
#         trees <- read.tree(paste0('c:/research/r/rogue-ms/data-raw/simulations/1/all.bs'))
#         profvis::profvis(Roguehalla(trees[1:5], dropset = 1))
# })
