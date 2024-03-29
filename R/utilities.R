#' @importFrom fastmatch fmatch
.NeverDrop <- function(neverDrop, labels) {
  if (missing(neverDrop) || is.null(neverDrop) || length(neverDrop) == 0) {
    return(character(0))
  } else {
    if (is.character(neverDrop)) {
      ndInt <- fmatch(neverDrop, labels)
      if (any(is.na(ndInt))) {
        warning("Can't keep taxa not in tree:\n  ",
                paste0(neverDrop[is.na(ndInt)], collapse = ", "))
      }
      neverDrop[!is.na(ndInt)]
    } else if (!is.numeric(neverDrop)) {
      stop("`neverDrop` must be of mode character or numeric")
    } else {
      if (any(neverDrop < 1)) {
        warning("Values of `neverDrop` < 1 ignored")
      }
      if (any(neverDrop > length(labels))) {
        warning("Values of `neverDrop` should correspond to leaves of tree")
      }
      labels[neverDrop[neverDrop > 1 & neverDrop <= length(labels)]]
    }
  }
}

.PrepareTrees <- function(trees) {
  # Edge lengths slow rearrangement without affecting calculations
  trees <- lapply(trees, function (tr) {tr$edge.length <- NULL; tr})
  trees <- lapply(trees, RenumberTips, trees[[1]])
  trees <- structure(lapply(trees, Preorder), class = "multiPhylo")
}
