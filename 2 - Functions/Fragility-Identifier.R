fragility_identifier <- function(x) {
  rle_x <- rle(x)
  # Identifying sequences of 1s with length >= 5
  idx <- which(rle_x$values == 1 & rle_x$lengths >= 5)
  if (length(idx) > 0) {
    end_positions <- cumsum(rle_x$lengths)[idx]
    start_positions <- end_positions - rle_x$lengths[idx] + 1
    return(list(start = start_positions, end = end_positions))
  } else {
    return(NULL)
  }
}
