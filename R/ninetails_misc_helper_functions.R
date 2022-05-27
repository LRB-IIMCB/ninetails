#' Winsorizes nanopore signal - removes high cliffs (extending above & below
#' signal range). Slightly affects signal extremes.
#'
#' Based on the A. Signorell's DescTools
#' https://github.com/AndriSignorell/DescTools/blob/master/R/DescTools.r
#'
#'
#' @param signal numeric vector. A vector corresponding to given ONT signal.
#'
#' @return winsorized signal (numeric vector).
#' @export
#'
#' @examples
#' \dontrun{
#'
#' winsorized_signal <- winsorize_signal(signal= sample(200:300))
#'
#' }
#'
winsorize_signal <- function(signal){

  #assertions

  if (missing(signal)) {
    stop("Signal vector is missing. Please provide a valid signal argument.", .call = FALSE)
  }

  assertthat::assert_that(assertive::is_numeric(signal),
                          msg=paste0("Signal vector must be numeric. Please provide a valid argument."))
  assertthat::assert_that(assertive::is_atomic(signal),
                          msg=paste0("Signal vector must be atomic. Please provide a valid argument."))
  assertthat::assert_that(assertive::is_atomic(signal),
                          msg=paste0("Signal vector must be atomic. Please provide a valid argument."))
  assertthat::assert_that(assertthat::noNA(signal),
                          msg="Signal vector must not contain missing values. Please provide a valid argument.")

  signal_q <- stats::quantile(x=signal, probs=c(0.005, 0.995), na.rm=TRUE, type=7)
  minimal_val <- signal_q[1L]
  maximal_val <- signal_q[2L]

  signal[signal<minimal_val] <- minimal_val
  signal[signal>maximal_val] <- maximal_val

  winsorized_signal <- as.integer(signal)

  return(winsorized_signal)

}


#' Counts overlapping segments created per given nanopore signal.
#'
#' @param signal numeric vector corresponding to the given ONT signal.
#' @param segment numeric [1]. Length of the chunk(s) to be created.
#' @param overlap numeric [1]. Length of the overlap between the chunks.
#'
#' @return number of overlapping chunks per given tail.
#' @export
#'
#' @examples
#'\dontrun{
#'
#' count_overlaps(signal = signal_vector,
#'                segment = 100
#'                overlap = 10)
#'}
#'
#'
count_overlaps <- function(signal, segment, overlap) {
  # avoiding 'no visible binding for global variable' error
  x <- NULL

  #Assertions
  if (missing(signal)) {
    stop("Signal vector is missing. Please provide a valid signal argument.", .call = FALSE)
  }

  if (missing(segment)) {
    stop("Segment number is missing. Please provide a valid segment argument.", .call = FALSE)
  }

  if (missing(overlap)) {
    stop("Overlap vector is missing. Please provide a valid overlap argument.", .call = FALSE)
  }

  assertthat::assert_that(assertive::is_numeric(signal),
                          msg=paste0("Signal vector must be numeric. Please provide a valid argument."))
  assertthat::assert_that(assertive::is_atomic(signal),
                          msg=paste0("Signal vector must be atomic. Please provide a valid argument."))
  assertthat::assert_that(overlap<segment,
                          msg="Segment value must be greater than overlap value. Please provide a valid argument.")

  checkmate::assert_integerish(segment,
                               tol = sqrt(.Machine$double.eps),
                               lower = 1,
                               upper = Inf,
                               any.missing = TRUE,
                               all.missing = TRUE,
                               len = NULL,
                               min.len = 1L,
                               max.len = 1,
                               unique = FALSE,
                               sorted = FALSE,
                               names = NULL,
                               typed.missing = FALSE,
                               null.ok = FALSE,
                               coerce = FALSE,
                               .var.name = checkmate::vname(x),
                               add = NULL
                               )

  checkmate::assert_integerish(overlap,
                               tol = sqrt(.Machine$double.eps),
                               lower = 0,
                               upper = Inf,
                               any.missing = TRUE,
                               all.missing = TRUE,
                               len = NULL,
                               min.len = 1L,
                               max.len = 1,
                               unique = FALSE,
                               sorted = FALSE,
                               names = NULL,
                               typed.missing = FALSE,
                               null.ok = FALSE,
                               coerce = FALSE,
                               .var.name = checkmate::vname(x),
                               add = NULL)





  starts <- seq(1, length(signal), by=segment-overlap)
  ends   <- starts + segment - 1

  counts_overlaps <- length(lapply(1:length(starts), function(i) signal[starts[i]:ends[i]]))

  return(counts_overlaps)
}


#' Loads keras model for multiclass signal prediction.
#'
#' @param keras_model_path either missing or character string. Full path of the .h5 file
#' with keras model used to predict signal classes. If function is called without this
#' argument(argument is missing) the default pretrained model will be loaded. Otherwise,
#' the dir with custom model shall be provided.
#'
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
    path_to_default_model <- system.file("extdata", "cnn_model", "model_grey_gasf100x100_32_64_updated.h5", package="ninetails")
    keras_model <- keras::load_model_hdf5(path_to_default_model)
  } else {
    keras_model <- keras::load_model_hdf5(keras_model_path)
  }
}



#' Checks if the provided directory contains fast5 files in the correct format.
#'
#' This function analyses the structure of the first fast5 file in the given
#' directory and checks whether it fulfills the analysis requirements (if the
#' file is multifast5, basecalled by Guppy basecaller and containing provided
#' basecall_group). Otherwise the function throws an error (with description).
#'
#' @param workspace character string. Full path of the directory to search the
#' basecalled fast5 files in. The Fast5 files have to be multi-fast5 file.
#'
#' @param basecall_group character string ["Basecall_1D_000"]. Name of the
#' level in the Fast5 file hierarchy from which the data should be extracted.
#'
#' @return outputs the text info with basic characteristics of the data.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' check_fast5_filetype <- function(workspace = '/path/to/guppy/workspace',
#'                                  basecalled_group = 'Basecall_1D_000')
#'
#' }
#'
#'
#'
# This lookup function is inspired by adnaniazi's explore-basecaller-and-fast5type.R from tailfindr
# https://github.com/adnaniazi/tailfindr/blob/master/R/explore-basecaller-and-fast5type.R

check_fast5_filetype <- function(workspace, basecall_group){

  #Assertions
  if (missing(workspace)) {
    stop("Directory with basecalled fast5s is missing. Please provide a valid workspace argument.", call. =FALSE)
  }

  if (missing(basecall_group)) {
    stop("Basecall group is missing. Please provide a valid basecall_group argument.", call. =FALSE)
  }

  assertthat::assert_that(assertive::is_character(workspace), msg = paste0("Path to fast5 files is not a character string. Please provide valid path to basecalled fast5 files."))


  #list fast5 files in given dir
  fast5_files_list <- list.files(path = workspace, pattern = "\\.fast5$", recursive = TRUE, full.names = TRUE)

  #count fast5 files
  num_fast5_files <- length(fast5_files_list)

  cat(paste0('[', as.character(Sys.time()), '] ','Found ', num_fast5_files, ' fast5 file(s) in provided directory.\n'))

  #closer look into the first file on the list
  selected_fast5_file <-fast5_files_list[1]
  selected_fast5_file_structure <- rhdf5::h5ls(file.path(selected_fast5_file), recursive = FALSE)

  selected_fast5_read <- selected_fast5_file_structure$name[1]

  cat(paste0('[', as.character(Sys.time()), '] ','Analyzing one of the given fast5 files to check', '\n','if the data are in required format... \n'))


  # checking whether fast5 file is single or multi
  is_multifast5 <- function(selected_fast5_file_structure){
    sum(grepl('read_', selected_fast5_file_structure$name)) > 0
  }

  assertthat::on_failure(is_multifast5) <- function(call, env) {
    paste0("The provided fast5 is single fast5 file. Please provide multifast5 file(s).")
  }

  assertthat::assert_that(is_multifast5(selected_fast5_file_structure))

  # check whether file is basecalled or not
  tryCatch(selected_basecall_group <- rhdf5::h5read(selected_fast5_file,paste0(selected_fast5_read,"/Analyses/", basecall_group)), error = function(e) { cat("The previewed fast5 file is not a basecalled one. Ninetails requires fast5 files basecalled by Guppy.") })

  if (exists('selected_basecall_group')) {
    # checking whether the fast5 file contains RNA ONT reads
    is_RNA <- function(selected_fast5_file, selected_fast5_read){
      read_context_tags <- rhdf5::h5readAttributes(selected_fast5_file,paste0(selected_fast5_read,"/context_tags"))
      read_context_tags$experiment_type == "rna"
    }

    assertthat::on_failure(is_RNA) <- function(call, env) {
      paste0("The provided fast5 does not contain RNA reads. Please provide multifast5 file(s) with RNA reads.")
    }
    assertthat::assert_that(is_RNA(selected_fast5_file, selected_fast5_read))


    # retrieve basecaller & basecalling model (read attributes)
    selected_basecall_group <- rhdf5::h5readAttributes(selected_fast5_file,paste0(selected_fast5_read,"/Analyses/", basecall_group))
    basecaller_used <- selected_basecall_group$name
    model_used <- selected_basecall_group$model_type

    # retrieve guppy basecaller version (read attributes)
    path_to_guppy_version <- rhdf5::h5readAttributes(selected_fast5_file,paste0(selected_fast5_read,"/tracking_id"))
    guppy_version <- path_to_guppy_version$guppy_version

    # close all handled instances (groups, attrs) of fast5 file
    rhdf5::h5closeAll()

    cat('  Previewed fast5 file parameters:\n')
    cat('    data type: RNA \n')
    cat('    fast5 file type: multifast5 \n')
    cat('    basecaller used:',basecaller_used,' \n')
    cat('    basecaller version:',guppy_version,' \n')
    cat('    basecalling model:',model_used,' \n',' \n')

  } else {
    # close all handled instances (groups, attrs) of fast5 file
    rhdf5::h5closeAll()
    stop(paste0('[', as.character(Sys.time()), '] ','Ninetails encountered an error. Please provide fast5 files basecalled by Guppy software. '))
  }

}


