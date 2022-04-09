#' Creates list of overlapping tail fragments extracted from
#' nanopolish output and fast5 files. Only fragments containing move==1
#' are included.
#'
#' @param tail_feature_list list object produced by create_tail_feature_list function.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#' @return A list of tail range chunks organized by the read ID
#' is returned. Only tails containing at least one move value equal to 1
#' are included. Always assign this returned list to a variable, otherwise
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

  # initial assertions

  if (missing(num_cores)) {
    stop("Number of declared cores is missing. Please provide a valid num_cores argument.", call. =FALSE)
  }

  if (missing(tail_feature_list)) {
    stop("List of features is missing. Please provide a valid tail_feature_list argument.", call. =FALSE)
  }

  assertthat::assert_that(assertive::is_numeric(num_cores), msg = paste("Declared core number must be numeric. Please provide a valid argument."))
  assertthat::assert_that(assertive::is_list(tail_feature_list),
                          msg = paste("Given tail_feature_list is not a list (class). Please provide valid file format."))

  # creating cluster for parallel computing
  doParallel::registerDoParallel(cores = num_cores)

  # first remove reads with only zero moved tails
  tail_feature_list <- Filter(function(x) sum(x$tail_moves) !=0, tail_feature_list)


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

    # work on subsets of reads in parallel
    tail_chunk_list <- c(tail_chunk_list, foreach::foreach(nam = names(tail_feature_list)[index_list[[indx]]])
                         %dopar% split_with_overlaps_moved(nam, tail_feature_list, segment = 100, overlap = 50))

    utils::setTxtProgressBar(pb, indx)

  }


  #flatten the list
  tail_chunk_list <- lapply(rapply(tail_chunk_list, enquote, how="unlist"), eval)

  close(pb)

  return(tail_chunk_list)

}
