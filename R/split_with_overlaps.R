# custom function to split vec - signal, move vec, base trace prob etc. into the list of chunks
# it adds NAs to fill the last chunk; nas are then turned into mean value
# based on https://stackoverflow.com/questions/8872376/split-vector-with-overlapping-samples-in-r

split_with_overlaps <- function(signal, segment, overlap) {
  starts <- seq(1, length(signal), by=segment-overlap)
  ends   <- starts + segment - 1

  split_signal <- lapply(1:length(starts), function(i) signal[starts[i]:ends[i]])

  # replace NAs with mean signal value
  result_split <- lapply(split_signal, function(n) replace(n, is.na(n), as.integer(mean(signal))))

  return(result_split)
}
