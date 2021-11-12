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
#'
#' @return
#' @export
#'
#' @examples 
#'\dontrun{
#'
#' check_fast5_filetype <- function(workspace = '/path/to/guppy/workspace',
#'                                  basecalled_group = 'Basecall_1D_000')
#'}
#'

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
  
  
  cat(paste0('Checking if given fast5 files are in required format.\n'))
  
  # checking whether fast5 file is single or multi
  is_multifast5 <- function(selected_fast5_file_structure){
    sum(grepl('read_', selected_fast5_file_structure$name)) > 0
  }
  
  assertthat::on_failure(is_multifast5) <- function(call, env) {
    paste0("The provided fast5 is single fast5 file. Please provide multifast5 file(s).")
  }
  
  assertthat::assert_that(is_multifast5(selected_fast5_file_structure))
  
  
  cat(paste0('Provided multifast5 file.\n'))
  
  # checking whether the selected basecall group is present within the selected file
  
  basecall_group_exists <- function(basecall_group, selected_fast5_file, selected_fast5_read){
    selected_basecall_group <- rhdf5::h5read(selected_fast5_file,paste0(selected_fast5_read,"/Analyses/", basecall_group))
    exists("selected_basecall_group")
  }
  
  assertthat::on_failure(basecall_group_exists) <- function(call, env) {
    paste0('Provided basecall_group is absent from the previewed fast5 file. Please provide correct basecall_group argument.')
  }
  
  assertthat::assert_that(basecall_group_exists(basecall_group, selected_fast5_file, selected_fast5_read))  
  cat(paste0('Checked fast5 file fulfills analysis requirements.\n')) 
  
  # close all handled instances (groups, attrs) of fast5 file
  rhdf5::h5closeAll()
  
  return(TRUE)
  
} # check_fast5_filetype

