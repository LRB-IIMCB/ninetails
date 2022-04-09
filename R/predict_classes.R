
#' Classification of matrices produced from read chunks with CNN.
#'
#' @param gasf_list [list] A list of gasf matrices organized by the read ID_index.
#'
#' @return a list of gasfs predictions based on used model.
#' @export
#'
#' @examples
#'\dontrun{
#'
#' predict_classes(gasf_list = gasf_list, keras_model = "/path/to/the_model.h5")
#'
#'}
predict_classes <- function(gasf_list){

  chunknames <- names(gasf_list)
  names(gasf_list)  <- NULL

  gasf_list <- simplify2array(gasf_list)
  gasf_list <- aperm(gasf_list, c(4,1,2,3))

  #normalise the data
  gasf_list <- gasf_list/255

  # Output info
  cat(paste('Classifying gramian angular summation fields...', '\n', sep=''))


  keras_model <- load_keras_model()

  #predict chunk class
  predicted_gasf_classes <- keras_model %>% predict(gasf_list) %>% keras::k_argmax()
  predicted_gasf_classes <- as.numeric(predicted_gasf_classes)

  predicted_list = list() # creating empty list for the extracted  data

  predicted_list[["chunkname"]] <- chunknames
  predicted_list[["prediction"]] <- predicted_gasf_classes

  # Output info
  cat(paste('Tensorflow finished predictions.', '\n', sep=''))


  return(predicted_list)
}
