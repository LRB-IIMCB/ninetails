#' Creates list of gramian angular summation matrices produced based on
#' list of splitted tails (tail chunks).
#'
#'
#' @param tail_chunk_list character string. Full path of the list object produced
#' by create_chunk_list function.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#'
#' @return A list of gasf matrices organized by the read ID_index
#' is returned. Always assign this returned list to a variable, otherwise
#' the long list will be printed to the console, which may crash your R session.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' create_gasf_list(feature_list = '/path/to/chunk_list', num_cores = 10)
#'
#'}
#'
#'

create_gasf_list <- function(tail_chunk_list, num_cores){

  # Assertions

  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }

  if (missing(tail_chunk_list)) {
    stop("List of tail chunks is missing. Please provide a valid chunk_list argument.", call. =FALSE)
  }

  assertthat::assert_that(assertive::is_numeric(num_cores), msg=paste("Declared core number must be numeric. Please provide a valid argument."))

  # avoiding 'no visible binding for global variable' error
  readname <- NULL

  # creating cluster for parallel computing
  doParallel::registerDoParallel(cores = num_cores)

  # this is list of indexes required for parallel computing; the main list is split for chunks
  index_list = split(1:length(names(tail_chunk_list)), ceiling(1:length(names(tail_chunk_list))/100))

  # header for progress bar
  cat(paste('Computing gramian angular summation fields...', '\n', sep=''))

  # progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(index_list), style = 3, width = 50, char = "=")

  #create empty list for extracted data
  gasf_list = list()

  # loop for parallel extraction
  for (indx in 1:length(index_list)){
    # use selected number of cores
    doParallel::registerDoParallel(cores = num_cores)

    # work on subsets of signals in parallel
    gasf_list <- c(gasf_list, foreach::foreach(nam = names(tail_chunk_list)[index_list[[indx]]]) %dopar% create_gasf(nam,tail_chunk_list))



    utils::setTxtProgressBar(pb, indx)

  }

  close(pb)

  #label each signal according to corresponding read name to avoid confusion
  gasf_names <- names(tail_chunk_list)
  names(gasf_list) <- gasf_names



  return(gasf_list)

}
