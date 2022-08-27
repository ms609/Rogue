#' Rogue
#'
#' "Rogue" implements approaches to identify rogue taxa in phylogenetic
#' analysis.
#' Rogues are wildcard leaves whose uncertain position, perhaps a result of
#' missing or conflicting data, reduces the resolution of consensus trees
#' \insertCite{Kearney2002}{Rogue}.
#' Consensus trees that omit rogue taxa can be more informative.
#'
#' "Rogue" allows the user to select a concept of "information" by which the
#' quality of consensus trees should be evaluated, and a heuristic approach
#' by which rogue taxa should be identified.
#'
#' Rogue detection using the phylogenetic and clustering information content
#' measures (\acronym{SPIC}, \acronym{SCIC}) \insertCite{SmithCons}{Rogue}
#' is implemented using a quick heuristic that drops
#' the least "stable" leaves one at a time,
#' using an _ad hoc_ definition of stability \insertCite{SmithCons}{Rogue};
#' and by a more exhaustive (and time-consuming) approach that considers
#' dropping all possible sets of up to _n_ leaves
#' \insertCite{Aberer2013}{Rogue}.
#'
#' The latter heuristic is implemented for the relative bipartition
#' "information" content and Pattengale's criterion
#' _via_ [RogueNaRok](https://rnr.h-its.org/about) \insertCite{Aberer2013}{Rogue}.
#'
#'
#' ## Citing "Rogue"
#'
#' If you find this package useful in your work, Please consider citing
#' Smith (2021).
#'
#' To cite the underlying methods, please cite
#' \insertCite{Aberer2013;textual}{Rogue} ("RogueNaRok")
#' or \insertCite{SmithCons;textual}{Rogue} (SPIC, SCIC), as appropriate.
#'
#'
#' @references
#' \insertAllCited{}
#'
#' @keywords internal
"_PACKAGE"

# Suppress "NOTE: Nothing imported from Rdpack":
#' @importFrom Rdpack reprompt
NULL
