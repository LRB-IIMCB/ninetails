#' Creates list of signal dataframes based on informations extracted from
#' nanopolish output and fast5 files.
#' 
#'
#' @param feature_list character string. Full path of the list object produced 
#' by create_feature_list function. 
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#'
#' @return A list of signal dataframes organized by the read ID 
#' is returned. Always assign this returned list to a variable, otherwise 
#' the long list will be printed to the console, which may crash your R session.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' create_signal_list(feature_list = '/path/to/feature_list',
#'                    num_cores = 10)
#' 
#'}
#'
#'

create_signal_list <- function(feature_list, num_cores){
  
  # Assertions
  
  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }
  
  if (missing(feature_list)) {
    stop("List of features is missing. Please provide a valid feature_list argument.", call. =FALSE)
  }
  
  assertthat::assert_that(assertive::is_numeric(num_cores), msg=paste("Declared core number must be numeric. Please provide a valid argument."))
  
  # CREATING CLUSTER FOR PARALLEL COMPUTING
  doParallel::registerDoParallel(cores = num_cores)
  
  # this is list of indexes required for parallel computing; the main list of reads is split for chunks 
  index_list = split(1:length(names(feature_list)), ceiling(1:length(names(feature_list))/100))
  
  #create empty list for extracted fast5 data
  signal_df_list = list()
  
  # loop for parallel extraction
  for (indx in 1:length(index_list)){
    # use selected number of cores
    doParallel::registerDoParallel(cores = num_cores)
    
    # work on subsets of reads in parallel
    signal_df_list <- c(signal_df_list, foreach::foreach(nam = names(feature_list)[index_list[[indx]]]) %dopar% create_signal_dataframe(nam,feature_list))
    
  }  
  
  #label each signal according to corresponding read name to avoid confusion
  squiggle_names <- names(feature_list)
  names(signal_df_list) <- squiggle_names
  
  return(signal_df_list)
  
}
