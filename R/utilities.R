.NeverDrop <- function (neverdrop, labels) {
  if (missing(neverdrop)) {
    return(integer(0))
  } else {
    if (is.character(neverdrop)) {
      ndInt <- match(neverdrop, labels)
      if (any(is.na(ndInt))) {
        warning("Can't keep taxa not in tree:\n  ",
                paste0(neverdrop[is.na(ndInt)], collapse = ', '))
      }
      ndInt[!is.na(ndInt)]
    } else if (!is.numeric(neverdrop)) {
      stop("`neverdrop` must be of mode character or numeric")
    } else {
      neverdrop
    }
  }
}
