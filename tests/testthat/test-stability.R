library('TreeTools', quietly = TRUE)

test_that("Errors handled", {
  expect_error(TipInstability(BalancedTree(4)))
  expect_error(TipInstability(as.phylo(0:3, 6), dev = 'error'))
  expect_error(TipInstability(as.phylo(0:3, 6), ave = 'error'))
})

test_that("Cophenetic() works", {
  Test <- function (tr) {
    tr <- Preorder(tr)
    tr$edge.length <- rep(1, nrow(tr$edge))
    tips <- seq_along(tr$tip.label)
    expect_equal(unname(ape::dist.nodes(tr)[tips, tips]), Cophenetic(tr))
    expect_equal(log(Cophenetic(tr)), Cophenetic(tr, log = TRUE))
  }
  Test(BalancedTree(4))
  Test(BalancedTree(6))
  Test(PectinateTree(7))
  Test(CollapseNode(BalancedTree(101), 104:111))
  Test(as.phylo(201, 1201))
})

test_that("ColByStability()", {
  expect_error(ColByStability(list(BalancedTree(7), BalancedTree(8))))
  trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')
  tipCol <- col2rgb(ColByStability(trees[3:6]))
  # plot(consensus(trees[3:6], p = 0.5), tip.col = ColByStability(trees[3:6]))
  expect_gt(tipCol[1, 'Rogue'], tipCol['red', 't1'])
  expect_equal(tipCol['t2'], tipCol['t1'])
})

