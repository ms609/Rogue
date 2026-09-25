> **Before starting work in this directory, read [`../AGENTS.md`](../AGENTS.md)**
> for multi-agent coordination rules, build/test infrastructure, GHA workflows,
> and worktree discipline. That file is the authoritative reference for all
> cross-package agent operations.

# Rogue - Agent Notes

## Package Overview

**Rogue** is an R package for identifying "rogue" (wildcard) taxa in
phylogenetic tree sets. Rogue taxa have uncertain positions that reduce consensus
tree resolution; removing them can increase information content. The package
provides information-theoretic detection methods (Smith 2022) and an interface to
the RogueNaRok C library (Aberer et al. 2013).

- **Language**: en-GB (British English throughout)
- Version, dependencies and system requirements: see `DESCRIPTION`.
  Changes are recorded in `NEWS.md`.

## Key Exported Functions

| Function | Purpose |
|----------|---------|
| `RogueTaxa()` | Main entry point — dispatches to method based on `info` param |
| `QuickRogue()` | Fast greedy heuristic (SPIC/SCIC via `TipInstability`) |
| `TipInstability()` | Per-leaf instability score (MAD/SD of graph geodesics) |
| `TipVolatility()` | Per-leaf volatility via phylogenetic info distance |
| `GraphGeodesic()` | Shortest-path distance matrix between leaves (wraps C) |
| `ColByStability()` | Colour vector for plotting leaf stability |
| `C_RogueNaRok()` | Direct interface to RogueNaRok C library (no input checks!) |
| `Cophenetic()` | Deprecated alias for `GraphGeodesic()` |

### Method dispatch in `RogueTaxa()`

- **rbic** → `.RogueNaRok()` (C library)
- **spic / scic** → `Roguehalla()` (exhaustive combinatorial, R)
- **fspic / fscic** → `QuickRogue()` (greedy heuristic, R + C)

## Architecture & Call Graph

```
RogueTaxa()
  ├─ .PrepareTrees()           # strip edge lengths, renumber tips, preorder
  ├─ .RogueNaRok()             # writes trees to tempfile, calls C, parses output
  ├─ Roguehalla()              # exhaustive: tests all dropset combinations
  └─ QuickRogue()              # greedy: iteratively drops least stable leaf
       └─ TipInstability()     # called once per iteration
            ├─ (fused path)    # TIP_INSTABILITY — one .Call returns scores
            └─ (fallback)      # per-tree GraphGeodesic() + matrixStats
```

**Performance-critical path**: `QuickRogue` → `TipInstability` →
`TIP_INSTABILITY` (`src/tip_instability.cpp`). One C++ call computes the
geodesics for every tree and reduces them to a per-leaf score, never
materialising the pairs × trees distance matrix as an R object. It requires
`log = TRUE` and the same number of edges in every tree (so polytomies can
disqualify a tree set); otherwise `TipInstability()` falls back to per-tree
`GraphGeodesic()` with `matrixStats` row statistics in R.

`QuickRogue()` is greedy and breaks ties with `which.max()`, so small changes
to instability numerics (median convention, MAD constant, pair ordering) can
change which taxa it drops. Keep the fused and fallback paths numerically
identical, and check vignette output after touching either.

## Source Layout

```
R/
  RogueTaxa.R        # RogueTaxa() dispatcher
  SPIC.R             # QuickRogue(), Roguehalla()
  stability.R        # TipInstability(), GraphGeodesic(), ColByStability(), TipVolatility()
  utilities.R        # .NeverDrop(), .PrepareTrees()
  zz_RogueNaRok.R    # .RogueNaRok(), C_RogueNaRok()

src/
  graph_geodesic.c   # geodesic worker, log lookup table, GRAPH_GEODESIC entry points
  tip_instability.cpp # TIP_INSTABILITY: fused geodesic + instability reduction
  geodesic.h         # declarations shared between the C and C++ units
  Rogue_init.c       # .Call registration
  Makevars           # Explicit C_SOURCES / CXX_SOURCES lists
  rnr/               # RogueNaRok C library (git submodule from ms609/RogueNaRok)

inst/example/150.bs  # Bootstrap tree file (150 taxa) used for benchmarking
```

## C/C++ Code Details

- `graph_geodesic_phylo()` (in `graph_geodesic.c`) is the core geodesic
  algorithm, adapted from `ape::dist.nodes()`; O(n²) in the number of nodes.
  Its inner loop (the `SETBOTH` macro) dominates profiles; this is fundamental
  to the algorithm and not easily improvable without a redesign.
- Log-transformed distances use the static lookup table `lg[]`, filled at load
  time and shared with the C++ unit through `geodesic.h`.
- `RogueNaRok` wraps the `rnr/` library.
- `.Call` entry points are registered in `Rogue_init.c` and called by symbol
  (`R_forceSymbols`), so a new entry point must be registered there before R
  can reach it.

## Workflow

- After modifying a function signature or roxygen documentation, always run
  `devtools::document()` followed by `devtools::check_man()`.
- After writing or editing documentation text, run
  `spelling::spell_check_package()`.

## Build Notes

- **Never use `devtools::load_all()`** for performance work — it compiles with
  `-O0`, which makes benchmarks meaningless.
- Use `pkgbuild::build()` or `R CMD INSTALL` for release-quality builds.
- When profiling with VTune, build with `-O2 -g -fno-omit-frame-pointer` and
  set `MAKEFLAGS="DLLFLAGS=-static-libgcc"` in Makevars.win.
- Clean stale `.o` files before switching between debug/release builds.
- `src/Makevars` lists all sources explicitly (no wildcard). When adding a
  source file, add it to `C_SOURCES` or `CXX_SOURCES` and register any new
  entry point in `Rogue_init.c`.

## Key Design Decisions

- **`.prepared` parameter**: `QuickRogue()` and `Roguehalla()` accept
  `.prepared = TRUE` so `RogueTaxa()` can call `.PrepareTrees()` once and
  pass the result through without redundant re-preparation.
- **Lower-triangle optimisation**: Distance matrices are symmetric, so
  instability is computed over the n(n-1)/2 unique leaf pairs only, then each
  pair's deviation is folded back onto both of its leaves.
- **`parallel` argument**: retained in `TipInstability()` and `QuickRogue()`
  for backwards compatibility only; the fused implementation is
  single-threaded and ignores it.

## Testing

- Test suite uses testthat.
- Example data: `inst/example/150.bs` (bootstrap trees, 150 taxa).
- Synthetic test trees generated via TreeTools (`BalancedTree`,
  `PectinateTree`, `AddTipEverywhere`).
- `vignettes/Bayesian.Rmd` downloads MrBayes output from `ms609/hyoliths`
  and falls back to synthetic trees when offline.
