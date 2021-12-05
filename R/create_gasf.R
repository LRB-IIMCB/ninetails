#' Converting ONT signal to gramian angular summation field.
#'
#' @param chunkname character string. Name of the given signal chunk within the
#' analyzed dataset.
#'
#' @param tail_chunk_list list object produced by create_tail_chunk_list function.
#'
#' @return an array (100,100,3) with values (RGB channels) representing ONT signal.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' create_gasf(chunkname = '1234-jhgfdr54-io643gju-ouy4378989', tail_chunk_list = tail_chunk_list)
#'
#'}
#'
#'
#'
create_gasf <- function(chunkname, tail_chunk_list){

  signal <- tail_chunk_list[[chunkname]]

  # rescale values so that all of them fall in the interval [-1, 1]:
  rescaled_signal <- (signal - max(signal) + (signal - min(signal))) / (max(signal) - min(signal))

  # calculate phi coefficient for interpolation to polar coordinates
  rescaled_signal <- acos(rescaled_signal)

  # create matrix by replicating vec
  rescaled_signal <- cbind(replicate(length(rescaled_signal), rescaled_signal))

  #calculate sum of phi
  rescaled_signal <- rescaled_signal + t(rescaled_signal)

  #calculate cosinus
  gasf_matrix <- round(cos(rescaled_signal), 4)

  #resize gasf to 100x100 matrix
  gasf_matrix <- resize_gasf(gasf_matrix, c(100,100))

  # assign HEX color values to gasf matrix
  gasf_matrix <- paintmap::color_matrix(gasf_matrix, colors=rainbow(100))

  #split color values to RGB channels
  gasf_matrix <- grDevices::col2rgb(gasf_matrix, alpha = FALSE)

  #reshape the data into new dimensions
  gasf_matrix <- array(t(gasf_matrix), c(100,100,3))


  return(gasf_matrix)
}
