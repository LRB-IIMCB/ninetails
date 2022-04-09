#' Loads keras model for multiclass signal prediction.
#'
#' @param keras_model_path either missing or character string. Full path of the .h5 file
#' with keras model used to predict signal classes. If function is called without this
#' argument(argument is missing) the default pretrained model will be loaded. Otherwise,
#' the dir with custom model shall be provided.
#'
#' @return
#' @export
#'
#' @examples
#'\dontrun{
#'
#' load_keras_model(keras_model_path = "/path/to/the/model/in/hdf5_format")
#' }
#'
load_keras_model <- function(keras_model_path){
  if (rlang::is_missing(keras_model_path)) {
    path_to_default_model <- system.file("extdata", "model_grey_gasf100x100_32_64.h5", package="ninetails")
    keras_model <- keras::load_model_hdf5(path_to_default_model)
  } else {
    keras_model <- keras::load_model_hdf5(keras_model_path)
  }
}
