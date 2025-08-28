test_that("RoguePlot basic functionality", {
  skip_if_not_installed("ape")
  library(ape)
  
  # Create test trees with rogue taxon
  trees <- list(
    read.tree(text = "(a, (b, (c, (rogue, (d, (e, f))))));"),
    read.tree(text = "(rogue, (a, (b, (c, (d, (e, f))))));"),
    read.tree(text = "((rogue, a), (b, (c, (d, (e, f)))));")
  )
  
  # Test basic function call
  result <- RoguePlot(trees, "rogue")
  
  expect_type(result, "list")
  expect_named(result, c("legendLabels", "frequencies", "consensus"))
  expect_type(result$legendLabels, "character")
  expect_type(result$frequencies, "double")
  expect_s3_class(result$consensus, "phylo")
})

test_that("RoguePlot error handling", {
  skip_if_not_installed("ape")
  library(ape)
  
  trees <- list(
    read.tree(text = "(a, (b, (c, (rogue, (d, (e, f))))));"),
    read.tree(text = "(rogue, (a, (b, (c, (d, (e, f))))));")
  )
  
  # Test missing arguments
  expect_error(RoguePlot(), "Both 'trees' and 'rogueTaxon' arguments are required")
  expect_error(RoguePlot(trees), "Both 'trees' and 'rogueTaxon' arguments are required")
  
  # Test invalid rogue taxon
  expect_error(RoguePlot(trees, "nonexistent"), "Rogue taxon 'nonexistent' not found in trees")
  
  # Test invalid input types
  expect_error(RoguePlot("not_a_list", "rogue"), "'trees' must be a list or multiPhylo object")
  
  # Test insufficient trees
  expect_error(RoguePlot(trees[1], "rogue"), "At least two trees are required")
})

test_that("RoguePlot with user example", {
  skip_if_not_installed("ape")
  library(ape)
  
  # User's exact example from issue
  trees <- list(
    read.tree(text = "(a, (b, (c, (rogue, (d, (e, f))))));"),
    read.tree(text = "(a, (b, (c, (rogue, (d, (e, f))))));"),
    read.tree(text = "(a, (b, (c, (rogue, (d, (e, f))))));"),
    read.tree(text = "(a, (b, (c, (rogue, (d, (e, f))))));"),
    read.tree(text = "(rogue, (a, (b, (c, (d, (e, f))))));"),
    read.tree(text = "((rogue, a), (b, (c, (d, (e, f)))));"),
    read.tree(text = "(a, (b, ((c, d), (rogue, (e, f)))));"),
    read.tree(text = "(a, (b, ((c, (rogue, d)), (e, f))));"),
    read.tree(text = "(a, (b, (c, (d, (rogue, (e, f)))))))")
  )
  
  # Test function with legend
  result <- RoguePlot(trees, "rogue", legend = "topleft", legend.inset = 0.02)
  
  expect_type(result, "list")
  expect_true(length(result$legendLabels) > 0)
  expect_true(all(result$frequencies >= 0))
  expect_true(all(result$frequencies <= 1))
})