#' Winsorizes signal, extracts poly(A) tail range and splits it for overlapping chunks.
#'
#' @param readname character string. Name of the given read within the
#' analyzed dataset.
#'
#' @param feature_list list object produced by create_feature_list function.
#'
#' @return a list of transformed poly(A) tail data: winsorized and splitted with overlaps.
#' @export
#'
#' @examples
#' \dontrun{
#' preprocess_tail(readname = 'abc123de-fg45-6789-0987-6543hijk2109', feature_list = feature_list)
#'
#'}
#'
preprocess_tail <- function(readname, feature_list){

  #assertion
  if (missing(readname)) {
    stop("Readname is missing. Please provide a valid readname argument.", call. =FALSE)
  }

  if (missing(feature_list)) {
    stop("List of reads' features is missing. Please provide a valid feature_list argument.", call. =FALSE)
  }

  assertthat::assert_that(assertive::is_list(feature_list), msg=paste("Given feature_list object is not a list. Please provide a valid argument."))
  assertthat::assert_that(readname %in% names(feature_list),msg = "Given feature_list does not contain provided readname. Please provide a valid readname argument.")

  #retrieve features
  polya_start_position <- feature_list[[readname]][[2]]
  transcript_start_position <- feature_list[[readname]][[3]]
  raw_signal <- feature_list[[readname]][[4]]

  #define polya end position
  polya_end_position <- transcript_start_position -1

  # winsorize the signal
  winsorized_signal <- winsorize_signal(raw_signal)

  # extract polya tail
  polya_tail <- winsorized_signal[polya_start_position:polya_end_position]

  #split to overlapping chunks
  splited_tail <- split_with_overlaps(polya_tail, 500, 200)

  return(splited_tail)


}
