#' Converting ONT signal to gramian angular summation field.
#'
#' @param chunkname character string. Name of the given signal chunk within the
#' analyzed dataset.
#'
#' @param tail_chunk_list list object produced by create_tail_chunk_list function.
#'
#' @return an array (100,100,1) with values (greyscale) representing ONT signal.
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

  #assertions
  if (missing(chunkname)) {
    stop("Chunkname is missing. Please provide a valid chunkname argument.", call. =FALSE)
  }

  if (missing(tail_chunk_list)) {
    stop("List of tail signal chunks is missing. Please provide a valid tail_chunk_list argument.", call. =FALSE)
  }

  assertthat::assert_that(assertive::is_character(chunkname),
                          msg = "Given chunkname is not a character string. Please provide a valid chunkname.")
  assertthat::assert_that(assertive::is_list(tail_chunk_list),
                          msg = "Given tail_chunk_list is not a list (class). Please provide valid file format.")


  tail_chunk <- tail_chunk_list[[chunkname]]

  # rescale values so that all of them fall in the interval [-1, 1]:
  tail_chunk <- (tail_chunk - max(tail_chunk) + (tail_chunk - min(tail_chunk))) / (max(tail_chunk) - min(tail_chunk))

  # calculate phi coefficient for interpolation to polar coordinates
  tail_chunk <- acos(tail_chunk)

  # create matrix by replicating vec
  tail_chunk <- cbind(replicate(length(tail_chunk), tail_chunk))

  #calculate sum of phi
  tail_chunk <- tail_chunk + t(tail_chunk)

  #calculate cosinus
  tail_chunk <- round(cos(tail_chunk), 4)

  # assign HEX color values to gasf matrix
  #tail_chunk <- paintmap::color_matrix(tail_chunk, colors=grDevices::rainbow(100))

  #split color values to RGB channels
  #tail_chunk <- grDevices::col2rgb(tail_chunk, alpha = FALSE)

  #reshape the data into new dimensions
  #tail_chunk <- array(t(tail_chunk), c(100,100,3))
  tail_chunk <- array(t(tail_chunk), c(100,100,1))


  return(tail_chunk)
}
