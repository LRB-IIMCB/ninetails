#' Creates vector of overlapping tail fragments extracted from
#' nanopolish output and fast5 files, in which move value = 1 occurs.
#'
#'
#' @param moves either a numeric vector of moves or character name
#' of such vector.
#'
#' @param signal either a corresponding numeric vector of signal values
#' or a character name of such vector.
#'
#' @param segment numeric [1]. Length of the chunk(s) to be created.
#'
#' @param overlap numeric [1]. Length of the overlap between the chunks.
#'
#' @return A vector of indices of moved chunks organized by their order in
#' an original list is returned.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' split_with_overlaps_moved(moves = vector_of_moves,
#'                           signal = vector_of_signal_corresponding_to_moves,
#'                           segment = 100
#'                           overlap = 10)
#'
#'}
#'
#'
split_with_overlaps_moved <- function(moves, signal, segment, overlap) {

  #assertions
  assertthat::assert_that(assertive::is_numeric(moves), msg=paste("Move vector must be numeric. Please provide a valid argument."))
  assertthat::assert_that(assertive::is_numeric(signal), msg=paste("Signal vector must be numeric. Please provide a valid argument."))
  assertthat::assert_that(assertive::is_atomic(signal), msg=paste("Signal vector must be atomic. Please provide a valid argument."))
  assertthat::assert_that(assertive::is_numeric(segment), msg=paste("Segment must be numeric. Please provide a valid argument."))
  assertthat::assert_that(assertive::is_numeric(overlap), msg=paste("Overlap must be numeric. Please provide a valid argument."))

  assertthat::assert_that(segment <= length(signal), msg="Overlap must not exceed the signal length. Please provide a valid argument.")
  assertthat::assert_that(segment>0, msg= "Segment value must be greater than 0. Please provide a valid argument.")

  assertthat::assert_that(overlap>=0, msg="Overlap value must be larger than or equal to 0. Please provide a valid argument.")
  assertthat::assert_that(overlap<segment, msg="Segment value must be greater than overlap value. Please provide a valid argument.")

  #prevent from running the function on zero-only vectors
  if (sum(moves) == 0){
    stop("The move vector does not contain the value of 1. Please provide a valid move argument (use create_tail_chunk_list_moved() function or prefilter data manually).", call. =FALSE)
  }

  #initial coordinates (for all chunks)
  start_coordinates_total <- seq(1, length(moves), by=segment-overlap)
  end_coordinates_total   <- start_coordinates_total + segment - 1

  #extract indices of "moved" chunks
  moved_chunks_indices <- lapply(1:length(start_coordinates_total), function(i) 1 %in% moves[start_coordinates_total[i]:end_coordinates_total[i]])
  moved_chunks_indices <- which(unlist(moved_chunks_indices))

  #coordinates of selected "moved" chunks
  start_coordinates_selected <- start_coordinates_total[moved_chunks_indices]
  end_coordinates_selected <- start_coordinates_selected + segment - 1

  #extract ONLY "moved" signal chunks
  extracted_moved_signals <- lapply(1:length(start_coordinates_selected), function(i) signal[start_coordinates_selected[i]:end_coordinates_selected[i]])
  # replace NAs with 3 most frequent values (randomly sampled)
  # if all values would be equal, so the breaks would not be unique

  most_freq_vals <- as.numeric(names(sort(table(signal),decreasing=TRUE)[1:3]))
  extracted_moved_signals <- lapply(extracted_moved_signals, function(n) replace(n, is.na(n), sample(most_freq_vals,sum(is.na(n)),TRUE)))

  #create selected chunk names from their indices
  names(extracted_moved_signals) <- moved_chunks_indices

  return(extracted_moved_signals)

}
