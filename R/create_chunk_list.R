#' Creates list of tails extracted based on informations from
#' nanopolish output and fast5 files.
#'
#'
#' @param feature_list list object produced by create_feature_list function.
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
#' create_chunk_list(feature_list = '/path/to/feature_list', num_cores = 10)
#'
#'}
#'
#'

create_chunk_list <- function(feature_list, num_cores){

  # Assertions

  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }

  if (missing(feature_list)) {
    stop("List of features is missing. Please provide a valid feature_list argument.", call. =FALSE)
  }

  assertthat::assert_that(assertive::is_numeric(num_cores), msg=paste("Declared core number must be numeric. Please provide a valid argument."))

  # creating cluster for parallel computing
  doParallel::registerDoParallel(cores = num_cores)

  # this is list of indexes required for parallel computing; the main list of reads is split for chunks
  index_list = split(1:length(names(feature_list)), ceiling(1:length(names(feature_list))/100))

  # header for progress bar
  cat(paste('Preprocessing poly(A) tails...', '\n', sep=''))

  # progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(index_list), style = 3, width = 50, char = "=")

  #create empty list for extracted data
  tail_chunk_list = list()

  # loop for parallel extraction
  for (indx in 1:length(index_list)){
    # use selected number of cores
    doParallel::registerDoParallel(cores = num_cores)

    # work on subsets of reads in parallel
    tail_chunk_list <- c(tail_chunk_list, foreach::foreach(nam = names(feature_list)[index_list[[indx]]]) %dopar% preprocess_tail(nam,feature_list))

    setTxtProgressBar(pb, indx)

  }

  close(pb)

  #label each signal according to corresponding read name to avoid confusion
  squiggle_names <- names(feature_list)
  names(tail_chunk_list) <- squiggle_names

  #naming chunks based on readnames & indices
  for (chunk in seq_along(tail_chunk_list)) {
    names(tail_chunk_list[[chunk]]) <- paste(names(tail_chunk_list)[chunk], seq_along(tail_chunk_list[[chunk]]), sep = "_")
  }

  list_structure <- rapply(tail_chunk_list, class)
  chunk_names <- names(list_structure)
  chunk_names <- gsub(".*\\.", "", chunk_names)

  #flatten the list
  tail_chunk_list <- Reduce(c, tail_chunk_list)

  names(tail_chunk_list) <- chunk_names


  return(tail_chunk_list)

}
