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

  assertthat::assert_that(assertive::is_character(workspace), msg = "Path to fast5 files is not a character string. Please provide valid path to basecalled fast5 files.")


  #list fast5 files in given dir
  fast5_files_list <- list.files(path = workspace, pattern = "\\.fast5$", recursive = TRUE, full.names = TRUE)

  #count fast5 files
  num_fast5_files <- length(fast5_files_list)

  cat(paste0('Found ', num_fast5_files, ' fast5 file(s) in provided directory.\n'))

  #closer look into the first file on the list
  selected_fast5_file <-fast5_files_list[1]
  selected_fast5_file_structure <- rhdf5::h5ls(file.path(selected_fast5_file), recursive = FALSE)

  selected_fast5_read <- selected_fast5_file_structure$name[1]

  cat(paste0('Analyzing one of the given fast5 files to check if the data are in required format... \n'))


  # checking whether fast5 file is single or multi
  is_multifast5 <- function(selected_fast5_file_structure){
    sum(grepl('read_', selected_fast5_file_structure$name)) > 0
  }

  assertthat::on_failure(is_multifast5) <- function(call, env) {
    paste0("The provided fast5 is single fast5 file. Please provide multifast5 file(s).")
  }

  assertthat::assert_that(is_multifast5(selected_fast5_file_structure))


  # checking whether the defined basecall group is present within the selected file
  basecall_group_exists <- function(basecall_group, selected_fast5_file, selected_fast5_read){
    selected_basecall_group <- rhdf5::h5read(selected_fast5_file,paste0(selected_fast5_read,"/Analyses/", basecall_group))
    exists("selected_basecall_group")
  }

  assertthat::on_failure(basecall_group_exists) <- function(call, env) {
  }

  assertthat::assert_that(basecall_group_exists(basecall_group, selected_fast5_file, selected_fast5_read))

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

  cat('Previewed fast5 file parameters:\n')
  cat('    data type: RNA \n')
  cat('    fast5 file type: multifast5 \n')
  cat('    basecaller used:',basecaller_used,' \n')
  cat('    basecaller version:',guppy_version,' \n')
  cat('    basecalling model:',model_used,' \n')

}
