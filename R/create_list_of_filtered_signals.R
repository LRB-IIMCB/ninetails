#' Creating list object with extracted filtered ONT signals.
#'
#' @param tail_chunk_list character string. Full path of the list object produced
#' by create_chunk_list function.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#' @return a list object with filtered names of signals containing peaks
#' and/or valleys.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' create_list_of_filtered_signals(tail_chunk_list='tail_chunk_list', num_cores=10)
#'}
#'
create_list_of_filtered_signals <- function(tail_chunk_list, num_cores, value){

  # Assertions
  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }

  if (missing(tail_chunk_list)) {
    stop("List of tail chunks is missing. Please provide a valid chunk_list argument.", call. =FALSE)
  }

  assertthat::assert_that(assertive::is_numeric(num_cores), msg=paste("Declared core number must be numeric. Please provide a valid argument."))

  # creating cluster for parallel computing
  doParallel::registerDoParallel(cores = num_cores)

  # this is list of indexes required for parallel computing; the main list is split for chunks
  index_list = split(1:length(names(tail_chunk_list)), ceiling(1:length(names(tail_chunk_list))/100))

  # header for progress bar
  cat(paste('Finding potentially modified chunks by signal threshold...', '\n', sep=''))

  # progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(index_list), style = 3, width = 50, char = "=")

  #create empty list for extracted data
  mod_list = list()

  # loop for parallel extraction
  for (indx in 1:length(index_list)){
    # use selected number of cores
    doParallel::registerDoParallel(cores = num_cores)

    # work on subsets of signals in parallel
    mod_list <- c(mod_list, foreach::foreach(nam = names(tail_chunk_list)[index_list[[indx]]]) %dopar% filter_by_threshold(nam, tail_chunk_list,lag = 50, threshold = 3, influence=0 ))

    setTxtProgressBar(pb, indx)

  }

  close(pb)

  #label each signal according to corresponding read name to avoid confusion
  mod_names <- names(tail_chunk_list)
  names(mod_list) <- mod_names

  #empirically tested on G's, 1 must be present more than 4 times in a row
  output_only_filtered <- mod_list[sapply(mod_list,function(vec) any(with(rle(vec==value), lengths > 4 & values)))]

  #output filtered read names instead of list of numbers
  filtered_names <- names(output_only_filtered)
  

  return(filtered_names)

}
