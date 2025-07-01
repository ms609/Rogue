# Rogue v2.1.6.9001 (2025-04-14)

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
