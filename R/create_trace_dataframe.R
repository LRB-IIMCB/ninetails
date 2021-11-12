#' Creates trace dataframe based on information retrieved from fast5 file.
#'
#' @param readname character string. Name of the given read within the
#' analyzed dataset. 
#'
#' @param feature_list character string. List of reads' features extracted
#' from nanopolish polya output, sequencing_summary and basecalled fast5
#' files. 
#'
#' @return A dataframe containing base probability values, positions, moves, traces, 
#' read segmentation (fctr: adapter, polya, transcript body). 
#' @export
#'
#' @examples
#'\dontrun{
#'
#' create_trace_dataframe(readname = 'abc123de-fg45-6789-0987-6543hijk2109',
#'                        feature_list = 'path/to/feature/list')
#' 
#'}


create_trace_dataframe <- function(readname, feature_list){
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
  raw_signal <- feature_list[[readname]][[7]]
  called_events <- feature_list[[readname]][[9]]
  stride <- feature_list[[readname]][[15]]
  trace <- feature_list[[readname]][[16]]
  
  # number of events expanded for whole signal vec (this is estimation of signal length, however keep in mind that decimal values are ignored) 
  number_of_events <- called_events * stride
  
  # convert trace to numeric
  trace_num <- as.numeric(trace) # values stored as coordinates with weird notation
  # create trace dataframe
  trace_df <- as.data.frame(matrix(trace_num, ncol = 8,  byrow = TRUE), stringsAsFactors = FALSE) # ncol taken from dims()
  # set colnames - there are no labels in fast5; adding labels makes table more visually pleasing
  colnames(trace_df) <- c("A_flip","C_flip","G_flip","T_flip","A_flop","C_flop","G_flop","T_flop")
  # add summarized base probabilities (I had to use "base" prefix, bc "T" is reserved for "TRUE" boolean arg in R)
  trace_df <- dplyr::mutate(trace_df, base_A = A_flip + A_flop, base_C = C_flip + C_flop, base_G = G_flip + G_flop, base_T = T_flip + T_flop)
  
  
  # Extract summarized base scores from trace_df
  baseA <- trace_df$base_A
  baseC <- trace_df$base_C
  baseG <- trace_df$base_G
  baseT <- trace_df$base_T
  
  # Create sample-wise vectors (to position traces along the signals)
  baseA_sample_wise <- c(rep(baseA, each=stride), rep(NA, length(raw_signal) - number_of_events))
  baseC_sample_wise <- c(rep(baseC, each=stride), rep(NA, length(raw_signal) - number_of_events))
  baseG_sample_wise <- c(rep(baseG, each=stride), rep(NA, length(raw_signal) - number_of_events))
  baseT_sample_wise <- c(rep(baseT, each=stride), rep(NA, length(raw_signal) - number_of_events))
  
  
  base_trace_df <- data.frame(baseA_sample_wise, baseC_sample_wise, baseG_sample_wise, baseT_sample_wise)
  
  #add position number; called events times stride -> length of this vec is equal to the called_events variable
  base_trace_df$position <- seq(from = 0, to = (nrow(base_trace_df)-1))
  
  return(base_trace_df)
  
}