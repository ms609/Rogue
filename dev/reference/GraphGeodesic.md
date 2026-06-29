# Graph Geodesic between leaves of unweighted tree

Graph Geodesic between leaves of unweighted tree

## Usage

``` r
GraphGeodesic(x, nTip = length(x$tip.label), log = FALSE, asMatrix = TRUE)

Cophenetic(x, nTip = length(x$tip.label), log = FALSE, asMatrix = TRUE)
```

## Arguments

- x:

  Object of class `phylo`.

- nTip:

  Integer specifying number of leaves.

- asMatrix:

  Logical specifying whether to coerce output to matrix format.

## Value

`GraphGeodesic()` returns an unnamed integer matrix describing the
number of edges between each pair of edges.

## Author

Martin R. Smith, modifying algorithm by Emmanuel Paradis in
[`ape::dist.nodes()`](https://rdrr.io/pkg/ape/man/cophenetic.phylo.html).

## Examples

``` r
GraphGeodesic(TreeTools::BalancedTree(5))
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    2    3    5    5
#> [2,]    2    0    3    5    5
#> [3,]    3    3    0    4    4
#> [4,]    5    5    4    0    2
#> [5,]    5    5    4    2    0
```
