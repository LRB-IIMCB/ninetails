#########################################################################
#   CREATE_SIGNAL_DATAFRAME 
#########################################################################

#' Creates signal dataframe based on retrieved features.
#' @param readname character string. Name of the given read within the
#' analyzed dataset. 
#'
#' @param feature_list character string. List of reads' features extracted
#' from nanopolish polya output, sequencing_summary and basecalled fast5
#' files. 
#'
#' @return A dataframe containing signal values, positions, moves, 
#' read segmentation (fctr: adapter, polya, transcript body). 
#' @export
#'
#' @examples
#'\dontrun{
#'
#' create_signal_dataframe(readname = 'abc123de-fg45-6789-0987-6543hijk2109',
#'                         feature_list = 'path/to/feature/list')
#' 
#'}
#'

create_signal_dataframe <- function(readname, feature_list){
  
  #assertion
  if (missing(readname)) {
    stop("Readname is missing. Please provide a valid readname argument.", call. =FALSE)
  }
  
  if (missing(feature_list)) {
    stop("List of reads' features is missing. Please provide a valid feature_list argument.", call. =FALSE)
  }
  
  assertthat::assert_that(assertive::is_list(feature_list), msg=paste("Given feature_list object is not a list. Please provide a valid argument."))
  assertthat::assert_that(readname %in% names(feature_list),msg = "Given feature_list does not contain provided readname. Please provide a valid readname argument.")
  
  # retrieve read features from feature list:
  called_events <- feature_list[[readname]][[9]]
  move <- feature_list[[readname]][[11]]
  polya_end_position <- feature_list[[readname]][[3]]
  polya_start_position <- feature_list[[readname]][[2]]
  raw_signal <- feature_list[[readname]][[7]]
  stride <- feature_list[[readname]][[15]]
  
  # number of events expanded for whole signal vec (this is estimation of signal length, however keep in mind that decimal values are ignored) 
  number_of_events <- called_events * stride
  # actual signal length (this allows to bind moves along signal, because it is with decimals)
  signal_length <- length(raw_signal) 
  
  #handling move data
  moves_sample_wise_vector <- c(rep(move, each=stride), rep(NA, signal_length - number_of_events))
  
  #handling event length
  #based on adnaniazi's tailfindr https://github.com/adnaniazi/tailfindr/blob/master/R/extract-read-data.R
  
  event_length_vector <- rep(0, length(move))
  count <- 0
  for (i in seq(from=called_events, to=1, by=-1)) {
    if (move[i] == 1) {
      event_length_vector[i] <- count + 1
      count <- 0
    } else {
      count <- count + 1
    }
  }
  
  #expand event_length_vector by stride
  event_length_vector <- event_length_vector * stride
  #expand event_length_vector by stride
  event_length_sample_wise_vector <- c(rep(event_length_vector, each=stride), rep(0, length(raw_signal) - number_of_events))
  # replace all subsequent non-NA values with 0 by rle function
  #based on StackOverflow wisdom, see:
  #https://stackoverflow.com/questions/38015358/replace-repeated-value-with-0-in-string
  event_length_sample_wise_vector <- event_length_sample_wise_vector *!duplicated(inverse.rle(within.list(rle(event_length_sample_wise_vector), values <-seq_along(values))))
  # reassign 0 as NAs
  event_length_sample_wise_vector[event_length_sample_wise_vector==0] <- NA
  
  #creating signal df
  signal_df <- data.frame(position=seq(1,signal_length,1), signal=raw_signal[1:signal_length], moves=moves_sample_wise_vector, event_length = event_length_sample_wise_vector) # this is signal converted to dframe
  
  # signal segmentation factor
  #adapter sequence
  signal_df$segment[signal_df$position < polya_end_position +1] <- "adapter"  
  #transcript sequence
  signal_df$segment[signal_df$position > polya_start_position] <- "transcript"
  #polya sequence
  signal_df$segment[signal_df$position>polya_start_position -1 & signal_df$position<polya_end_position +1] <- "poly(A)"
  
  
  signal_df$segment <- as.factor(signal_df$segment)
  
  return(signal_df)
  
}