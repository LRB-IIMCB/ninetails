#' Creates tail range dataframe based on retrieved features.
#' @param readname character string. Name of the given read within the
#' analyzed dataset. 
#'
#' @param feature_list character string. List of reads' features extracted
#' from nanopolish polya output, sequencing_summary and basecalled fast5
#' files. 
#'
#' @return A dataframe containing signal values, positions, moves, 
#' read segmentation (fctr: adapter, polya, transcript body) within
#' the tail range defined as poly(A) tail +/- 2 times stride. 
#' @export
#'
#' @examples
#'\dontrun{
#'
#' create_tail_dataframe(readname = 'abc123de-fg45-6789-0987-6543hijk2109',
#'                       feature_list = 'path/to/feature/list')
#' 
#'}
#'

create_tail_dataframe <- function(readname, feature_list){
  
  #assertion
  if (missing(readname)) {
    stop("Readname is missing. Please provide a valid readname argument.", call. =FALSE)
  }
  
  if (missing(feature_list)) {
    stop("List of reads' features is missing. Please provide a valid feature_list argument.", call. =FALSE)
  }
  
  assertthat::assert_that(assertive::is_list(feature_list), msg=paste("Given feature_list object is not a list. Please provide a valid argument."))
  assertthat::assert_that(readname %in% names(feature_list),msg = "Given feature_list does not contain provided readname. Please provide a valid readname argument.")
  
  #retrieve polya coordinates
  polya_end_position <- feature_list[[readname]][[3]]
  polya_start_position <- feature_list[[readname]][[2]]
  
  #construct signal dataframe
  signal_df <- create_signal_dataframe(readname, feature_list)
  
  #construct trace dataframe
  trace_df <- create_trace_dataframe(readname, feature_list)
  
  #combine dataframes
  merged_df <- dplyr::full_join(signal_df, trace_df, by="position")
  
  #subset signal dataframe to tail region
  #plus/minus 10 positions down- and upstream from tail
  
  tail_range_df <- subset(merged_df, position>polya_start_position -11 & position<polya_end_position +11)
  
  
  return(tail_range_df)
  
}