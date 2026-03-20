# Rogue - Agent Notes

## Package Overview

**Rogue** (v2.1.7) is an R package for identifying "rogue" (wildcard) taxa in
phylogenetic tree sets. Rogue taxa have uncertain positions that reduce consensus
tree resolution; removing them can increase information content. The package
provides information-theoretic detection methods (Smith 2022) and an interface to
the RogueNaRok C library (Aberer et al. 2013).

- **Language**: en-GB (British English throughout)
- **License**: GPL (>= 3)
- **System requirements**: C99
- **Author/maintainer**: Martin R. Smith

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
            ├─ (batch path)    # LOG_GRAPH_GEODESIC_MULTI — single .Call for all trees
            └─ (per-tree path) # GraphGeodesic() loop — fallback for heterogeneous trees
```

**Performance-critical path**: `QuickRogue` → `TipInstability` → C geodesic +
Rfast row statistics. The batch C path (`LOG_GRAPH_GEODESIC_MULTI`) computes
lower-triangle distances for all trees in one call, reusing a single interim
buffer. It requires `log = TRUE` and uniform tree dimensions; otherwise falls
back to per-tree `GraphGeodesic()`.

## Source Layout

```
R/
  RogueTaxa.R        # RogueTaxa() dispatcher
  SPIC.R             # QuickRogue(), Roguehalla()
  stability.R        # TipInstability(), GraphGeodesic(), ColByStability(), TipVolatility()
  utilities.R        # .NeverDrop(), .PrepareTrees()
  zz_RogueNaRok.R   # .RogueNaRok(), C_RogueNaRok()
  Rogue-package.R    # Package-level docs

src/
  graph_geodesic.c   # GRAPH_GEODESIC, LOG_GRAPH_GEODESIC, LOG_GRAPH_GEODESIC_MULTI
  Rogue_init.c       # .Call registration (4 C functions)
  Makevars           # Explicit SOURCES/OBJECTS list
  rnr/               # RogueNaRok C library (git submodule from ms609/RogueNaRok)

tests/testthat/
  test-RogueTaxa.R, test-stability.R, test-spic.R, test-RogueNaRok.R, test-utilities.R
  testdata/          # Fixture data
inst/example/
  150.bs             # Bootstrap tree file (150 taxa) used for benchmarking
```

## C Code Details

`graph_geodesic.c` contains four registered `.Call` functions:

1. **`GRAPH_GEODESIC`** (5 args) — integer distance matrix, all nodes × all nodes
2. **`LOG_GRAPH_GEODESIC`** (5 args) — log-transformed doubles, tip × tip only;
   uses a static lookup table (up to 32768)
3. **`LOG_GRAPH_GEODESIC_MULTI`** (6 args) — batch version of the log path;
   accepts concatenated edge arrays for all trees, returns lower-triangle only
   (n_pairs × n_trees), reuses single interim buffer across trees
4. **`RogueNaRok`** (10 args) — wraps the rnr/ library

The core algorithm (`graph_geodesic_phylo`) is O(n²) in the number of nodes.
VTune profiling shows the inner-loop `SETBOTH` macro (~line 56) accounts for
~62% of Rogue.dll time — this is fundamental to the algorithm and not easily
improvable without a redesign.

## Dependencies

**Key imports**: ape (phylo I/O), TreeTools (tree manipulation), TreeDist
(information measures), Rfast (fast row stats: `rowMads`, `rowMedians`,
`rowVars`, `rowmeans`), matrixStats (supplementary), fastmatch (fast matching),
cli (progress bars).

**Suggests**: testthat, knitr, rmarkdown, spelling, PlotTools.

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
- The `src/Makevars` file lists all C sources explicitly (no wildcard).
  When adding a new `.c` file, update both SOURCES and register in
  `Rogue_init.c`.

## Key Design Decisions

- **`.prepared` parameter**: `QuickRogue()` and `Roguehalla()` accept
  `.prepared = TRUE` so `RogueTaxa()` can call `.PrepareTrees()` once and
  pass the result through without redundant re-preparation.
- **Auto-parallel heuristic**: `TipInstability()` automatically enables OpenMP
  in Rfast operations when `nrow(dists_lt) > 1000L`, avoiding overhead on small
  problems.
- **Lower-triangle optimisation**: Distance matrices are symmetric, so
  `TipInstability()` operates on only the n(n-1)/2 unique pairs, then
  reconstructs the full symmetric deviation matrix for `rowmeans()`.
- **Batch C auto-fallback**: `LOG_GRAPH_GEODESIC_MULTI` requires all trees to
  have the same number of edges; if trees differ (e.g. polytomies),
  `TipInstability()` falls back to per-tree `GraphGeodesic()`.

## Testing

- Test suite uses testthat with snapshot testing (`_snaps/` directory).
- Example data: `inst/example/150.bs` (bootstrap trees, 150 taxa).
- Synthetic test trees generated via TreeTools (`BalancedTree`,
  `PectinateTree`, `AddTipEverywhere`).
