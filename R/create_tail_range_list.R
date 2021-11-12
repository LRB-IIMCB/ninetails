#' Creates list of tail range dataframes based on informations extracted from
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
#' @return A list of tail range dataframes organized by the read ID 
#' is returned. Always assign this returned list to a variable, otherwise 
#' the long list will be printed to the console, which may crash your R session.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' create_tail_range_list(feature_list = '/path/to/feature_list',
#'                        num_cores = 10)
#' 
#'}
#'
#'

create_tail_range_list <- function(feature_list, num_cores){
  
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
  
  # header for progress bar
  cat(paste('Creating list of tail ranges...', '\n', sep=''))
  
  # progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(index_list), style = 3, width = 50, char = "=")
  
  #create empty list for extracted fast5 data
  tail_range_list = list()
  
  # loop for parallel extraction
  for (indx in 1:length(index_list)){
    # use selected number of cores
    doParallel::registerDoParallel(cores = num_cores)
    
    # work on subsets of reads in parallel
    tail_range_list <- c(tail_range_list, foreach::foreach(nam = names(feature_list)[index_list[[indx]]]) %dopar% create_tail_dataframe(nam,feature_list))
    
    setTxtProgressBar(pb, indx)
    
  }  
  
  close(pb)
  
  #label each signal according to corresponding read name to avoid confusion
  squiggle_names <- names(feature_list)
  names(tail_range_list) <- squiggle_names
  
  return(tail_range_list)
  
}