.NeverDrop <- function (neverdrop, labels) {
  if (missing(neverdrop)) {
    return(character(0))
  } else {
    if (is.character(neverdrop)) {
      ndInt <- match(neverdrop, labels)
      if (any(is.na(ndInt))) {
        warning("Can't keep taxa not in tree:\n  ",
                paste0(neverdrop[is.na(ndInt)], collapse = ', '))
      }
      neverdrop[!is.na(ndInt)]
    } else if (!is.numeric(neverdrop)) {
      stop("`neverdrop` must be of mode character or numeric")
    } else {
      if (any(neverdrop < 1)) {
        warning("Values of `neverdrop` < 1 ignored")
      }
      if (any(neverdrop > length(labels))) {
        warning("Values of `neverdrop` should correspond to leaves of tree")
      }
      labels[neverdrop[neverdrop > 1 & neverdrop <= length(labels)]]
    }
  }
}
