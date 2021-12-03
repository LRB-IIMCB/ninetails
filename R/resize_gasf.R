#' Resizes gramian angular summation field.
#'
#' @param gasf - an array, gramian angular summation matrix produced by
#' create_gasf() function.
#'
#' @param new_dims numeric vector. Dimensions of new matrix.
#'
#' @return resized gasf (according to defined dims)
#' @export
#'
#' @examples
#' \dontrun{
#' resize_gasf <- function(gasf=gasf, new_dims=c(100,100,3))
#' }

#https://stackoverflow.com/questions/11123152/function-for-resizing-matrices-in-r

resize_gasf <- function(gasf, new_dims=dim(gasf)){

  # input gasf
  obj_dims <- dim(gasf)
  obj <- list(x= 1:obj_dims[1], y=1:obj_dims[2], z= gasf)

  # output gasf
  rescaled_gasf <- matrix(NA, nrow=new_dims[1], ncol=new_dims[2])
  new_dims <- dim(rescaled_gasf)

  rescale <- function(x, new_range=range(x)){
    x_range <- range(x)
    matrix_factor <- (new_range[2]-new_range[1])/(x_range[2]-x_range[1])
    new_range[1]+(x-x_range[1])*matrix_factor
  }

  # rescaling
  new_coordinates <- as.matrix(expand.grid(seq_len(new_dims[1]), seq_len(new_dims[2])))
  locations <- new_coordinates
  locations[,1] = rescale(new_coordinates[,1], c(1, obj_dims[1]))
  locations[,2] = rescale(new_coordinates[,2], c(1, obj_dims[2]))

  # interpolation
  rescaled_gasf[new_coordinates] <- fields::interp.surface(obj, locations)
  rescaled_matrix <- round(rescaled_gasf, 4)

  return(rescaled_matrix)
}

