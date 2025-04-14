library("TreeTools", quietly = TRUE)

test_that("TipInstability() error handling", {
  expect_error(TipInstability(as.phylo(0:3, 6), dev = "error"))
  expect_error(TipInstability(as.phylo(0:3, 6), ave = "error"))
})

test_that("ColByStability() error handling", {
  expect_equal(ColByStability(NULL), character(0))
})

test_that("TipInstability() null output", {
  expect_equal(TipInstability(BalancedTree(4)),
               setNames(rep(0, 4), paste0("t", 1:4)))
})

test_that("GraphGeodesic() works", {
  Test <- function(tr) {
    tr <- Preorder(tr)
    tr$edge.length <- rep(1, nrow(tr$edge))
    tips <- seq_along(tr$tip.label)
    expect_equal(unname(ape::dist.nodes(tr)[tips, tips]), GraphGeodesic(tr))
    expect_equal(log(GraphGeodesic(tr)), GraphGeodesic(tr, log = TRUE))
  }
  Test(BalancedTree(4))
  Test(BalancedTree(6))
  Test(PectinateTree(7))
  Test(CollapseNode(BalancedTree(101), 104:111))
  Test(as.phylo(201, 1201))
  
  expect_equal(GraphGeodesic(BalancedTree(4), log = FALSE),
               GraphGeodesic(BalancedTree(4), log = 1))
})

test_that("ColByStability()", {
  expect_error(ColByStability(list(BalancedTree(7), BalancedTree(8))))
  trees <- AddTipEverywhere(BalancedTree(8), "Rogue")
  tipCol <- col2rgb(ColByStability(trees[3:6]))
  # plot(consensus(trees[3:6], p = 0.5), tip.col = ColByStability(trees[3:6]))
  expect_gt(tipCol[1, "Rogue"], tipCol["red", "t1"])
  expect_equal(tipCol["t2"], tipCol["t1"])
})

test_that("ColByStability() - stable trees", {
  trees <- c(BalancedTree(4), BalancedTree(4), BalancedTree(4))
  cons <- Consensus(trees)
  tipCols <- ColByStability(trees)[cons$tip.label]
  expect_equal(tipCols, setNames(hcl.colors(131, "inferno")[rep(1, 4)],
                                 cons[["tip.label"]]))
})

