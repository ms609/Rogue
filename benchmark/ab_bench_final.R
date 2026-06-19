.libPaths(c("C:/Users/pjjg18/GitHub/Rogue/.ab_dev",
            "C:/Users/pjjg18/GitHub/Rogue/.ab_ref",
            .libPaths()))
library(RogueRef)
library(RogueDev)
library(TreeTools)
library(bench)

trees150 <- ape::read.tree(system.file("example/150.bs", package="RogueRef",
                           lib.loc="C:/Users/pjjg18/GitHub/Rogue/.ab_ref"))
trees_med <- trees150[1:100]
trees_small <- AddTipEverywhere(BalancedTree(8), "Rogue")

# Warmup (includes OpenMP thread pool init)
invisible(RogueRef::TipInstability(trees_small))
invisible(RogueDev::TipInstability(trees_small))
invisible(RogueRef::QuickRogue(trees_small, info="phylogenetic"))
invisible(RogueDev::QuickRogue(trees_small, info="phylogenetic"))

cat("=== Canary: GraphGeodesic (single tree, should be neutral) ===\n")
b <- bench::mark(
  ref = RogueRef::GraphGeodesic(trees_small[[1]], log=TRUE),
  dev = RogueDev::GraphGeodesic(trees_small[[1]], log=TRUE),
  min_iterations=50, check=FALSE)
print(b[, c("expression","min","median","n_itr")])

cat("\n=== TipInstability (100 trees x 150 tips) ===\n")
b <- bench::mark(
  ref = RogueRef::TipInstability(trees_med, log=TRUE, average="median", deviation="mad"),
  dev = RogueDev::TipInstability(trees_med, log=TRUE, average="median", deviation="mad"),
  min_iterations=5, check=FALSE)
print(b[, c("expression","min","median","n_itr","mem_alloc")])

cat("\n=== QuickRogue (100 trees x 150 tips) ===\n")
b <- bench::mark(
  ref = RogueRef::QuickRogue(trees_med, info="phylogenetic"),
  dev = RogueDev::QuickRogue(trees_med, info="phylogenetic"),
  min_iterations=3, check=FALSE)
print(b[, c("expression","min","median","n_itr","mem_alloc")])

cat("\n=== RogueTaxa fspic (100 trees x 150 tips) ===\n")
b <- bench::mark(
  ref = RogueRef::RogueTaxa(trees_med, info="fspic"),
  dev = RogueDev::RogueTaxa(trees_med, info="fspic"),
  min_iterations=3, check=FALSE)
print(b[, c("expression","min","median","n_itr","mem_alloc")])

# Verify correctness
cat("\n=== Correctness check ===\n")
r_ref <- RogueRef::QuickRogue(trees_med, info="phylogenetic")
r_dev <- RogueDev::QuickRogue(trees_med, info="phylogenetic")
cat("Results match:", identical(r_ref, r_dev), "\n")

ti_ref <- RogueRef::TipInstability(trees_med, log=TRUE, average="median", deviation="mad")
ti_dev <- RogueDev::TipInstability(trees_med, log=TRUE, average="median", deviation="mad")
cat("TipInstability match:", all.equal(ti_ref, ti_dev), "\n")
