# VTune driver: exercise the C hot paths in TipInstability
library(Rogue, lib.loc = ".vtune-lib")
library(TreeTools)

trees150 <- ape::read.tree(system.file("example/150.bs", package = "Rogue",
                                        lib.loc = ".vtune-lib"))
trees_med <- trees150[1:100]

# Warmup
invisible(TipInstability(trees_med, log = TRUE, average = "median",
                          deviation = "mad"))

# Repeat enough times for meaningful VTune sampling (~30s target)
cat("Starting VTune workload...\n")
for (i in seq_len(50)) {
  TipInstability(trees_med, log = TRUE, average = "median",
                  deviation = "mad", checkTips = FALSE)
}
cat("Done.\n")
