test_that("Cophenetic() works", {
  library(TreeTools)
  Test <- function (tr) {
    tr <- Preorder(tr)
    tr$edge.length <- rep(1, nrow(tr$edge))
    tips <- seq_along(tr$tip.label)
    expect_equal(unname(ape::dist.nodes(tr)[tips, tips]), Cophenetic(tr))
  }
  Test(BalancedTree(4))
  Test(BalancedTree(6))
  Test(PectinateTree(7))
  Test(CollapseNode(BalancedTree(101), 104:111))
})

test_that("ColByStability()", {
  expect_error(ColByStability(list(BalancedTree(7), BalancedTree(8))))
  trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')
  tipCol <- col2rgb(ColByStability(trees[3:6]))
  expect_lt(tipCol[1, 'Rogue'], tipCol['red', 't1'])
  expect_equal(tipCol['t2'], tipCol['t1'])
})

