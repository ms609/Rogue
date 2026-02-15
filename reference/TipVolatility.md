# Detect rogue taxa using phylogenetic information distance

Calculate the volatility of each tip: namely, the impact on the mean
phylogenetic information distance (Smith 2020) between trees when that
tip is removed. Effective when the number of trees is small.

## Usage

``` r
TipVolatility(trees)
```

## Arguments

- trees:

  List of trees to analyse.

## Value

`TipVolatility()` returns a named vector listing the volatility index
calculated for each leaf. Higher values indicate more volatile leaves.

## References

Smith MR (2020). “Information theoretic Generalized Robinson–Foulds
metrics for comparing phylogenetic trees.” *Bioinformatics*, **36**(20),
5007–5013.
[doi:10.1093/bioinformatics/btaa614](https://doi.org/10.1093/bioinformatics/btaa614)
.

## See also

Other tip instability functions:
[`TipInstability()`](https://ms609.github.io/Rogue/reference/TipInstability.md)

## Examples

``` r
library("TreeTools", quietly = TRUE)

# Generate some trees with two rogue taxa
trees <- AddTipEverywhere(BalancedTree(8), "Rogue")
trees[] <- lapply(trees, AddTip, "Rogue", "Rogue2")

# Calculate tip volatility
sb <- TipVolatility(trees)

# Use volatility to colour leaves in consensus tree
sbNorm <- 1 + (99 * (sb - min(sb)) / (max(sb - min(sb))))
col <- hcl.colors(128, "inferno")[sbNorm]
plot(consensus(trees), tip.color = col)

# Add a legend for the colour scale used
PlotTools::SpectrumLegend(
  "bottomleft", bty = "n", # Suppress box
  inset = -0.02,           # Avoid overlap
  title = "Volatility",
  legend = signif(seq(max(sb), min(sb), length.out = 4), 3),
  palette = hcl.colors(128, "inferno")
)


# Plot consensus after removing highly volatile taxa
plot(ConsensusWithout(trees, names(sb[sb == max(sb)])))
```
