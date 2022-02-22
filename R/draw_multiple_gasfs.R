#' Drawing multiple gramian angular summation matrices produced based on
#' list of splitted tails (tail chunks).
#'
#'
#' @param gasf_list character string. Full path of the list object produced
#' by create_gasf_list function.
#'
#' @param num_cores numeric [1]. Number of physical cores to use in processing
#' the data. Do not exceed 1 less than the number of cores at your disposal.
#'
#'
#' @return A set of plotted gasfs.
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' draw_multiple_gasfs(gasf_list = '/path/to/gasf_list', num_cores = 10)
#'
#'}
#'
#'

draw_multiple_gasfs <- function(gasf_list, num_cores){

  # creating cluster for parallel computing
  doParallel::registerDoParallel(cores = num_cores)

  # this is list of indexes required for parallel computing; the main list is split for chunks
  index_list = split(1:length(names(gasf_list)), ceiling(1:length(names(gasf_list))/100))

  # header for progress bar
  cat(paste('Drawing graphs...', '\n', sep=''))

  # progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(index_list), style = 3, width = 50, char = "=")

  #create empty list for extracted data
  plot_list = list()

  # loop for parallel extraction
  for (indx in 1:length(index_list)){
    # use selected number of cores
    doParallel::registerDoParallel(cores = num_cores)

    # work on subsets of signals in parallel
    plot_list <- c(plot_list, foreach::foreach(nam = names(gasf_list)[index_list[[indx]]]) %dopar% draw_gasf(nam,gasf_list))



    utils::setTxtProgressBar(pb, indx)

  }

  close(pb)

  #label each signal according to corresponding read name to avoid confusion
  gasf_names <- names(gasf_list)
  names(plot_list) <- gasf_names



  #return(plot_list)


}
