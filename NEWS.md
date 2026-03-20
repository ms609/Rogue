# Rogue v2.2.0 (2026-03-20)

## Performance

- `TipInstability()` operates on the lower triangle of distance matrices only,
  halving the work for row statistics (`rowMads`, `rowVars`, `rowMedians`).

- New batch C function `LOG_GRAPH_GEODESIC_MULTI` computes log-geodesic
  distances for all trees in a single `.Call()`, reusing one interim buffer
  and returning lower-triangle entries directly.

- Cache-friendly extraction loop in `graph_geodesic.c` (stride-1 access
  instead of stride-`all_nodes`).

- Auto-enable OpenMP parallelism in Rfast operations (`rowMads`, `rowVars`)
  when the distance matrix exceeds 1 000 rows.

- Use `Rfast::rowMedians()` in place of `matrixStats::rowMedians()`.

- `RogueTaxa()` now calls `.PrepareTrees()` once and passes `.prepared = TRUE`
  to `QuickRogue()` / `Roguehalla()`, avoiding redundant tree preparation.

- `QuickRogue()` precomputes information upper bounds for all iterations rather
  than recomputing each step.


# Rogue v2.1.7 (2025-07-01)

- Improve tip instability calculation in identical tree sets
  ([#29](https://github.com/ms609/Rogue/issues/29)).
  
- Improve variable protection.


# Rogue v2.1.6 (2023-11-29)

- Legend annotations in documentation.
- Disable parallel evaluation by default in `TipInstability()`,
  adding `parallel` parameter to allow user to override.
- Use format string in REprintf().


# Rogue v2.1.5 (2023-03-20)

- Call C functions using symbols, not strings.


# Rogue v2.1.4 (2023-01-16)

- C2X compliant function prototypes.

- Remove unused `sprintf()` calls.


# Rogue v2.1.3 (2022-09-26)

- `ColByStability()` gains `pal` argument to allow specification of custom
  palettes.


# Rogue v2.1.2 (2022-08-16)

- Faster rogue detection when edge lengths provided, per report by Joe Keating.

- Don't list `neverDrop` in `QuickRogue(fullSeq = TRUE)`.


# Rogue v2.1.1 (2022-07-20)

- Handle `ColByStability(trees = NULL)`.


# Rogue v2.1.0 (2022-01-13)

- Early termination of `QuickRogue()` when no further improvement possible.

- `Cophenetic()` renamed to the more accurate `GraphGeodesic()`.

- Calculate information content of consensus trees with p > 0.5 (#15).

- Improve support for `multiPhylo` objects.

- New vignette detailing rogue detection with Bayesian tree samples.


# Rogue v2.0.0 (2021-09-13)

- Information theoretic rogue detection (per Smith, 2022).


# Rogue v1.0.0 (2021-06-28)

 - R interface to RogueNaRok.
