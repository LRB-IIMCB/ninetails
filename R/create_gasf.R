#' Converting ONT signal to gramian angular summation field.
#'
#' @param tail_chunk character string. Name of the given signal chunk within the
#' analyzed dataset.
#'
#' @return an array (100,100,1) with values representing ONT signal.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' create_gasf(tail_chunk = "1234-jhgfdr54-io643gju-ouy4378989")
#'
#'}
#'
create_gasf <- function(tail_chunk){

  #assertions
  if (missing(tail_chunk)) {
    stop("Tail_chunk is missing. Please provide a valid tail_chunk argument.", call. =FALSE)
  }

  # rescale values so that all of them fall in the interval [-1, 1]:
  tail_chunk <- (tail_chunk-max(tail_chunk)+(tail_chunk-min(tail_chunk)))/(max(tail_chunk)-min(tail_chunk))

  # calculate phi coefficient for interpolation to polar coordinates
  tail_chunk <- acos(tail_chunk)

  # create matrix by replicating vec
  tail_chunk <- cbind(replicate(length(tail_chunk), tail_chunk))

  #calculate sum of phi
  tail_chunk <- tail_chunk + t(tail_chunk)

  #calculate cosinus
  tail_chunk <- cos(tail_chunk)

  #reshape the data into new dimensions
  tail_chunk <- array(t(tail_chunk), c(100,100,1))

  # rescale values so that all of them fall in the interval [0, 1]:
  tail_chunk <- round((tail_chunk-min(tail_chunk))/(max(tail_chunk)-min(tail_chunk)), 4)


  return(tail_chunk)
}
