#' Creates list of gramian angular summation matrices produced based on
#' list of splitted tails (tail chunks).
#'
#' @param tail_chunk_list character string. Full path of the list object produced
#' by create_chunk_list function.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
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
    stop("List of tail chunks is missing. Please provide a valid tail_chunk_list argument.", call. =FALSE)
  }

  assertthat::assert_that(assertive::is_numeric(num_cores), msg=paste("Declared core number must be numeric. Please provide a valid argument."))


  #register cores for parallelization
  doParallel::registerDoParallel(cores = num_cores)

  #retrieve chunknames
  chunknames <- gsub(".*?\\.","",names(rapply(tail_chunk_list, function(x) head(x, 1))))

  #create empty list for the data
  gasf_list = list()

  #set progressbar
  cat(paste('Computing gramian angular summation fields...', '\n', sep=''))
  pb <- utils::txtProgressBar(min = 0, max = length(tail_chunk_list), style = 3, width = 50, char = "=")

  #loop through the nested list
  for (read in seq_along(tail_chunk_list)){
    for (chunk in seq_along(tail_chunk_list[[read]])){
      gasf_list <- c(gasf_list, foreach::foreach(chunk) %dopar% create_gasf(tail_chunk_list[[read]][[chunk]]))

      utils::setTxtProgressBar(pb, read)
    }
  }

  close(pb)

  #restore names in the list
  names(gasf_list) <- chunknames

  return(gasf_list)

}



