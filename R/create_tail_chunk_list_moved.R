#' Creates list of overlapping tail fragments extracted from
#' nanopolish output and fast5 files. Only fragments containing move==1
#' are included.
#'
#'
#' @param tail_feature_list list object produced by create_tail_feature_list function.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#'
#' @return A list of tail range chunks organized by the read ID
#' is returned. Always assign this returned list to a variable, otherwise
#' the long list will be printed to the console, which may crash your R session.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' create_tail_chunk_list_moved(tail_feature_list = '/path/to/tail_feature_list', num_cores = 10)
#'
#'}
#'
#'


create_tail_chunk_list_moved <- function(tail_feature_list, num_cores){

  # Assertions

  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }

  if (missing(tail_feature_list)) {
    stop("List of features is missing. Please provide a valid tail_feature_list argument.", call. =FALSE)
  }

  assertthat::assert_that(assertive::is_numeric(num_cores), msg=paste("Declared core number must be numeric. Please provide a valid argument."))

  # creating cluster for parallel computing
  doParallel::registerDoParallel(cores = num_cores)

  # this is list of indexes required for parallel computing; the main list of reads is split for chunks
  index_list = split(1:length(names(tail_feature_list)), ceiling(1:length(names(tail_feature_list))/100))

  # header for progress bar
  cat(paste('Filtering chunks based on move values...', '\n', sep=''))

  # progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(index_list), style = 3, width = 50, char = "=")

  #create empty list for extracted data
  tail_chunk_list = list()

  # loop for parallel extraction
  for (indx in 1:length(index_list)){
    # use selected number of cores
    doParallel::registerDoParallel(cores = num_cores)

    # work on subsets of reads in parallel
    tail_chunk_list <- c(tail_chunk_list, foreach::foreach(nam = names(tail_feature_list)[index_list[[indx]]])
                         %dopar% split_with_overlaps_moved(moves = tail_feature_list[[nam]][[5]],
                                                           signal = tail_feature_list[[nam]][[4]],
                                                           segment = 100,overlap = 50))

    utils::setTxtProgressBar(pb, indx)

  }


  # label each signal according to corresponding read name to avoid confusion
  names(tail_chunk_list) <- names(tail_feature_list)


  # naming chunks based on readnames & indices
  chunk_names <- paste0(rep(names(tail_chunk_list), lengths(tail_chunk_list)), '_', unlist(lapply(tail_chunk_list, names)))


  #flatten the list
  tail_chunk_list <- Reduce(c, tail_chunk_list)

  # rename flattened list
  names(tail_chunk_list) <- chunk_names


  close(pb)


  return(tail_chunk_list)

}
